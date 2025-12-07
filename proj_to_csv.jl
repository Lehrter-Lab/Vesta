using DataFrames, CSV, Statistics, JSON3
using NCDatasets
using ArchGDAL
using ThreadsX, Base.Threads
using FilePathsBase

function compute_local_data(ncfile::String; timesteps_per_chunk::Int=10, grid_size::Float64=2000.0)
    println("DEBUG: Opening NetCDF file $ncfile"); flush(stdout)
    ds = NCDataset(ncfile, "r")
    N_records = length(ds["pid"])
    pid_all = ds["pid"][:]

    # -------------------------
    # Determine number of particles and timesteps
    # Assumption: pid is sequential from 1 to max(pid) and repeats each timestep
    # -------------------------
    wrap_idx = findfirst(diff(pid_all) .< 0)
    N_particles = wrap_idx
    t_steps = N_records / N_particles
    println("DEBUG: Detected $N_particles particles per timestep over $t_steps timesteps"); flush(stdout)

    # -------------------------
    # Determine chunk size based on timesteps_per_chunk
    # -------------------------
    chunk_size = N_particles * timesteps_per_chunk
    total_chunks = ceil(Int, N_records / chunk_size)
    println("DEBUG: Processing in $total_chunks chunks of up to $chunk_size records each"); flush(stdout)

    # -------------------------
    # Load coordinates
    # -------------------------
    x_all = ds["lon"][:]
    y_all = ds["lat"][:]
    time_all = ds["time"][:]

    min_x, max_x = extrema(x_all)
    min_y, max_y = extrema(y_all)
    edges_x = collect(min_x:grid_size:max_x)
    edges_y = collect(min_y:grid_size:max_y)
    n_x, n_y = length(edges_x)-1, length(edges_y)-1

    # Initialize buffers
    dt_sum_cell        = zeros(Float64, n_x, n_y)
    n_particles_cell   = zeros(Int,     n_x, n_y)

    nthreads = Threads.nthreads()
    master_dt = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    master_np = [zeros(Int,     n_x, n_y) for _ in 1:nthreads]

    # -------------------------
    # Process chunks
    # -------------------------
    for i in 0:(total_chunks-1)
        start_idx = i*chunk_size + 1
        stop_idx  = min((i+1)*chunk_size, N_records)
        println("DEBUG: Chunk $i/$total_chunks  range = $start_idx:$stop_idx"); flush(stdout)

        x_chunk    = x_all[start_idx:stop_idx]
        y_chunk    = y_all[start_idx:stop_idx]
        time_chunk = time_all[start_idx:stop_idx]

        n_chunk_records = stop_idx - start_idx + 1
        n_chunk_timesteps = ceil(Int, n_chunk_records / N_particles)

        # -------------------------
        # Process each particle instance sequentially
        # -------------------------
        @threads for p in 1:N_particles
            tid = threadid()
            local_dt = master_dt[tid]
            local_np = master_np[tid]

            # Track a particle’s current visit
            prev_bin_x = 0
            prev_bin_y = 0
            accum_dt = 0.0
            
            for t in 2:n_chunk_timesteps
                idx_prev = (t-2)*N_particles + p
                idx_curr = (t-1)*N_particles + p
                if idx_curr > n_chunk_records
                    break
                end
            
                dt_val = time_chunk[idx_curr] - time_chunk[idx_prev]
                if dt_val <= 0
                    continue
                end
            
                x_val = x_chunk[idx_curr]
                y_val = y_chunk[idx_curr]
                x_bin = searchsortedlast(edges_x, x_val) - 1
                y_bin = searchsortedlast(edges_y, y_val) - 1
            
                in_bounds = (1 ≤ x_bin ≤ n_x && 1 ≤ y_bin ≤ n_y)
            
                # If still in same cell: accumulate dt
                if in_bounds && x_bin == prev_bin_x && y_bin == prev_bin_y
                    accum_dt += dt_val
            
                # If entering a new valid cell:
                elseif in_bounds
                    # Finalize previous visit if one existed
                    if prev_bin_x != 0
                        @inbounds begin
                            local_dt[prev_bin_x, prev_bin_y] += accum_dt
                            local_np[prev_bin_x, prev_bin_y] += 1
                        end
                    end
                    # Start new visit
                    prev_bin_x = x_bin
                    prev_bin_y = y_bin
                    accum_dt = dt_val
            
                # If leaving the grid:
                else
                    if prev_bin_x != 0
                        @inbounds begin
                            local_dt[prev_bin_x, prev_bin_y] += accum_dt
                            local_np[prev_bin_x, prev_bin_y] += 1
                        end
                    end
                    prev_bin_x = 0
                    prev_bin_y = 0
                    accum_dt = 0.0
                end
            end
            
            # Finalize last visit at end of chunk
            if prev_bin_x != 0
                @inbounds begin
                    local_dt[prev_bin_x, prev_bin_y] += accum_dt
                    local_np[prev_bin_x, prev_bin_y] += 1
                end
            end
        end
    end

    close(ds)

    # -------------------------
    # Merge threads
    # -------------------------
    println("DEBUG: Merging thread results..."); flush(stdout)
    for tid in 1:nthreads
        dt_sum_cell        .+= master_dt[tid]
        n_particles_cell   .+= master_np[tid]
    end

    # -------------------------
    # Build DataFrame
    # -------------------------
    println("DEBUG: Total particle counts in grid: ", sum(n_particles_cell))
    println("DEBUG: Max dt_sum in grid: ", maximum(dt_sum_cell))
        
    rows = Vector{NamedTuple}(undef, n_x * n_y)
    k = 1
    for xi in 1:n_x, yi in 1:n_y
        dt_val = dt_sum_cell[xi, yi]
        np_val = n_particles_cell[xi, yi]
        mean_val = np_val > 0 ? dt_val / np_val : NaN
        total_val = np_val > 0 ? dt_val : NaN
    
        rows[k] = (x_bin = xi, y_bin = yi,
                   dt_sum = dt_val,
                   n_particles = np_val,
                   mean_exp_time = mean_val,
                   total_exp_time = total_val)
        k += 1
    end

    df = DataFrame(rows)

    csv_path = replace(ncfile, ".nc" => ".csv")
    CSV.write(csv_path, df)

    meta_path = replace(ncfile, ".nc" => ".meta.json")
    meta = Dict("grid_size" => grid_size,
                "bounds" => (min_x, max_x, min_y, max_y),
                "n_x" => n_x,
                "n_y" => n_y)
    JSON3.write(meta_path, meta)

    return df
end

# -------------------------
# Main function
# -------------------------
function main()
    println("DEBUG: Entering main()"); flush(stdout)

    input_line = try readline(stdin) catch _ "" end
    args = split(strip(input_line))
    println("DEBUG: Parsed args = $(args)"); flush(stdout)
  
    defaults = Dict("ncfile" => "particle_enhanced_proj.nc",
                "grid_size" => "2500.0",
                "timesteps_per_chunk" => "10")

    input_dict = Dict{String,String}()
    for arg in args
        if occursin('=', arg)
            k, v = split(arg, "=", limit=2)
            input_dict[strip(k)] = strip(v)
        end
    end

    params = merge(defaults, input_dict)
    ncfile = params["ncfile"]
    grid_size = parse(Float64, params["grid_size"])
    timesteps_per_chunk = parse(Int, params["timesteps_per_chunk"])
    println("DEBUG: Running with ncfile=$ncfile grid=$grid_size timesteps_per_chunk=$timesteps_per_chunk"); flush(stdout)

    csv_path  = replace(ncfile, ".nc" => ".csv")
    meta_path = replace(ncfile, ".nc" => ".meta.json")

    println("DEBUG: Running compute_local_data..."); flush(stdout)
    compute_local_data(ncfile; grid_size=grid_size, timesteps_per_chunk=timesteps_per_chunk)
    
    println("DEBUG: All done."); flush(stdout)
end

# Call main
main()

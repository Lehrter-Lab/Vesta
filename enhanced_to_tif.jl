using DataFrames, CSV, Statistics, JSON3
using ArchGDAL, GeometryBasics, GeoDataFrames, Proj, NCDatasets, GeoTables
using CairoMakie, Colors, ThreadsX, Base.Threads
using FilePathsBase

# Compute and Save Tabular Results (CSV + Metadata)
function compute_local_data(ncfile::String;
                            grid_size::Float64=2000.0,
                            target_crs::String="EPSG:5070",
                            chunk_size::Int=1_000_000)

    ds = NCDataset(ncfile, "r")
    N  = length(ds["pid"])

    # -------------------------
    # Bounding Box Checkpoint
    # -------------------------
    bbox_tmpfile = replace(ncfile, ".nc" => ".bbox.tmp")

    if isfile(bbox_tmpfile)
        println("Found bounding box checkpoint: $bbox_tmpfile")
        bbox_data        = readlines(bbox_tmpfile)
        min_lon, max_lon = parse.(Float64, split(bbox_data[1]))
        min_lat, max_lat = parse.(Float64, split(bbox_data[2]))
        println("Loaded bounding box from checkpoint:")
        println("  lon: [$min_lon, $max_lon], lat: [$min_lat, $max_lat]")
    else
        println("Pass 1: Computing domain bounding box...")
        min_lon, max_lon = Inf, -Inf
        min_lat, max_lat = Inf, -Inf

        for start in 1:chunk_size:N
            stop = min(start+chunk_size-1, N)
            lons = @views ds["lon"][start:stop]
            lats = @views ds["lat"][start:stop]

            good = .!(isnan.(lons) .| isnan.(lats) .| (abs.(lons) .> 1e3) .| (abs.(lats) .> 1e3))
            if any(good)
                min_lon = min(min_lon, minimum(lons[good]))
                max_lon = max(max_lon, maximum(lons[good]))
                min_lat = min(min_lat, minimum(lats[good]))
                max_lat = max(max_lat, maximum(lats[good]))
            end
        end

        if !isfinite(min_lon) || !isfinite(min_lat)
            close(ds)
            error("No valid lon/lat values found in $ncfile")
        end

        open(bbox_tmpfile, "w") do io
            println(io, "$min_lon $max_lon")
            println(io, "$min_lat $max_lat")
            flush(io)
        end
        println("Saved bounding box checkpoint to $bbox_tmpfile")
    end

    # -------------------------
    # Projection Setup
    # -------------------------
    trans_fwd = Proj.Transformation("EPSG:4326", target_crs; always_xy=true)
    trans_inv = Proj.Transformation(target_crs, "EPSG:4326"; always_xy=true)
    println("Created transformations successfully.")

    # Project Bounding Box
    corners = [(min_lon, min_lat), (max_lon, min_lat), (min_lon, max_lat), (max_lon, max_lat)]
    proj_corners = trans_fwd.(first.(corners), last.(corners))
    if isa(proj_corners, Tuple) && length(proj_corners) == 2 &&
       isa(proj_corners[1], AbstractArray) && isa(proj_corners[2], AbstractArray)
        xs, ys = proj_corners
    else
        xs = first.(proj_corners)
        ys = last.(proj_corners)
    end
    min_x, max_x = extrema(xs)
    min_y, max_y = extrema(ys)

    # Build Equal-Area Grid
    edges_x  = min_x:grid_size:max_x
    edges_y  = min_y:grid_size:max_y
    n_x, n_y = length(edges_x)-1, length(edges_y)-1

    # Initialize Buffers
    dt_sum_cell        = zeros(Float64, n_x, n_y)
    time_weighted_cell = zeros(Float64, n_x, n_y)
    n_particles_cell   = zeros(Int, n_x, n_y)

    nthreads = Threads.nthreads()
    thread_dt = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    thread_tw = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    thread_np = [zeros(Int,     n_x, n_y) for _ in 1:nthreads]

    # -------------------------
    # Progress checkpoint setup
    # -------------------------
    total_chunks = ceil(Int, N / chunk_size)
    chunk_indices = collect(1:chunk_size:N)
    progress_file = replace(ncfile, ".nc" => ".progress.tmp")

    last_chunk_done = 0
    if isfile(progress_file)
        last_chunk_done = parse(Int, readlines(progress_file)[1])
        println("Resuming from chunk $last_chunk_done")
    end

    chunks_per_10pct = max(1, ceil(Int, total_chunks / 10))

    # -------------------------
    # Chunk processing loop
    # -------------------------
    for i in (last_chunk_done + 1):total_chunks
        start_idx = chunk_indices[i]
        stop_idx  = min(start_idx+chunk_size-1, N)

        println("Processing chunk $i / $total_chunks")

        @views begin
            pid_chunk  = ds["pid"][start_idx:stop_idx]
            lon_chunk  = ds["lon"][start_idx:stop_idx]
            lat_chunk  = ds["lat"][start_idx:stop_idx]
            time_chunk = ds["time"][start_idx:stop_idx]
        end

        # Transform coordinates
        lon_vec = lon_chunk[:]
        lat_vec = lat_chunk[:]
        
        # Apply transformation using broadcasting
        x_vec, y_vec = trans_fwd.(lon_vec, lat_vec)
        
        # Reshape back to original chunk shape
        x_chunk = reshape(x_vec, size(lon_chunk))
        y_chunk = reshape(y_vec, size(lon_chunk))

        lon_chunk = nothing
        lat_chunk = nothing

        # Sort and group by pid
        sorted_idx  = sortperm(pid_chunk)
        pid_sorted  = pid_chunk[sorted_idx]
        x_sorted    = x_chunk[sorted_idx]
        y_sorted    = y_chunk[sorted_idx]
        time_sorted = time_chunk[sorted_idx]

        pid_chunk  = nothing
        x_chunk    = nothing
        y_chunk    = nothing
        time_chunk = nothing

        breaks = vcat(1, findall(diff(pid_sorted) .!= 0) .+ 1, length(pid_sorted)+1)

        @threads for g in 1:(length(breaks)-1)
            tid = threadid()
            local_dt = thread_dt[tid]
            local_tw = thread_tw[tid]
            local_np = thread_np[tid]

            lo = breaks[g]
            hi = breaks[g+1]-1

            if hi > lo
                for idx in (lo+1):hi
                    dt_val = time_sorted[idx] - time_sorted[idx-1]
                    if dt_val > 0
                        x_val = x_sorted[idx]
                        y_val = y_sorted[idx]
                        x_bin = clamp(Int(floor((x_val - edges_x[1]) / grid_size)) + 1, 1, n_x)
                        y_bin = clamp(Int(floor((y_val - edges_y[1]) / grid_size)) + 1, 1, n_y)
                        @inbounds begin
                            local_dt[x_bin, y_bin] += dt_val
                            local_tw[x_bin, y_bin] += dt_val * time_sorted[idx]
                            local_np[x_bin, y_bin] += 1
                        end
                    end
                end
            end
        end

        pid_sorted  = nothing
        x_sorted    = nothing
        y_sorted    = nothing
        time_sorted = nothing

        # -------------------------
        # Save progress checkpoint every ~10% or last chunk
        # -------------------------
        if i % chunks_per_10pct == 0 || i == total_chunks
            open(progress_file, "w") do io
                println(io, i)
                flush(io)
            end
            println("Saved progress checkpoint at chunk $i / $total_chunks")
        end
    end

    close(ds)

    # -------------------------
    # Merge threads
    # -------------------------
    for tid in 1:nthreads
        dt_sum_cell        .+= thread_dt[tid]
        time_weighted_cell .+= thread_tw[tid]
        n_particles_cell   .+= thread_np[tid]
    end

    # -------------------------
    # Save final CSV + metadata
    # -------------------------
    n_rows = count(!iszero, n_particles_cell)
    rows = Vector{NamedTuple}(undef, n_rows)
    k = 1
    for xi in 1:n_x, yi in 1:n_y
        if n_particles_cell[xi, yi] > 0
            rows[k] = (x_bin = xi, y_bin = yi,
                       dt_sum = dt_sum_cell[xi, yi],
                       time_weighted = time_weighted_cell[xi, yi],
                       n_particles = n_particles_cell[xi, yi])
            k += 1
        end
    end

    df = DataFrame(rows)
    df.mean_exp_time  = df.dt_sum ./ df.n_particles
    df.mean_water_age = df.time_weighted ./ df.dt_sum

    csv_path = replace(ncfile, ".nc" => ".csv")
    if isfile(csv_path)
        println("CSV already exists. Skipping recomputation: $csv_path")
    else
        CSV.write(csv_path, df)
        println("Saved results to $csv_path")
    end

    meta_path = replace(ncfile, ".nc" => ".meta.json")
    meta = Dict("grid_size" => grid_size,
                "target_crs" => target_crs,
                "bounds" => (min_x, max_x, min_y, max_y),
                "n_x" => n_x,
                "n_y" => n_y)
    JSON3.write(meta_path, meta)
    println("Saved compact metadata to $meta_path")

    # -------------------------
    # Remove progress checkpoint after successful completion
    # -------------------------
    try
        if isfile(progress_file)
            rm(progress_file; force=true)
            println("Deleted progress checkpoint: $progress_file")
        end
    catch err
        @warn "Could not delete progress checkpoint" exception=(err, catch_backtrace())
    end

    # Delete bbox checkpoint if present
    try
        if isfile(bbox_tmpfile)
            rm(bbox_tmpfile; force=true)
            println("Deleted bounding box checkpoint: $bbox_tmpfile")
        end
    catch err
        @warn "Could not delete bounding box checkpoint" exception=(err, catch_backtrace())
    end

    return df
end

function export_geospatial(csv_path::String, meta_path::String; fmt::String="GTiff")
    df = CSV.read(csv_path, DataFrame)
    meta = JSON3.read(Base.read(meta_path, String))

    grid_size    = meta["grid_size"]
    (min_x, max_x, min_y, max_y) = meta["bounds"]
    n_x, n_y     = meta["n_x"], meta["n_y"]
    target_crs   = meta["target_crs"]
    geotransform = [min_x, grid_size, 0.0, max_y, 0.0, -grid_size]
    numeric_cols = filter(c -> eltype(df[!, c]) <: Real && !(c in [:x_bin, :y_bin]), names(df))

    for col in numeric_cols
        output_path = replace(csv_path, ".csv" => "_$(col).tif")
        driver = ArchGDAL.getdriver(fmt)
        ArchGDAL.create(driver;
            filename=output_path,
            width=n_x,
            height=n_y,
            nbands=1,
            dtype=Float32) do dataset
                srs = ArchGDAL.importEPSG(parse(Int, split(target_crs, ":")[2]))
                ArchGDAL.setproj!(dataset, ArchGDAL.toWKT(srs))
                ArchGDAL.setgeotransform!(dataset, geotransform)

                band = ArchGDAL.getband(dataset, 1)
                data = fill(NaN32, n_y, n_x)
                for r in eachrow(df)
                    xi, yi = r.x_bin, r.y_bin
                    if 1 <= xi <= n_x && 1 <= yi <= n_y
                        data[n_y - yi + 1, xi] = Float32(r[col])
                    end
                end
                ArchGDAL.write!(band, data)
        end
        println("Exported single-band GeoTIFF → $output_path")
    end

    return [replace(csv_path, ".csv" => "_$(col).tif") for col in numeric_cols]
end

# Main
function main(; resume=false)
    # Example STDIN input:
    # echo "ncfile=particle_enhanced.nc grid_size=5000.0 crs_proj=EPSG:5070 chunk_size=10000000" | julia script.jl

    # Allowed keys and defaults
    valid_keys = Set(["ncfile", "grid_size", "crs_proj", "chunk_size"])
    defaults = Dict("ncfile"     => "particle_enhanced.nc",
                    "grid_size"  => "5000.0",
                    "crs_proj"   => "EPSG:5070",
                    "chunk_size" => "10000000")

    # Read optional line from STDIN
    input_line = try readline(stdin) catch; "" end
    args = split(strip(input_line))

    # Parse name=value pairs into a Dict
    input_dict = Dict{String,String}()
    for arg in args
        if occursin('=', arg)
            k, v = split(arg, "=", limit=2)
            k = strip(k); v = strip(v)

            if k ∈ valid_keys
                input_dict[k] = v
            else
                @warn "Unrecognized argument name: $k (ignored)"
            end
        elseif !isempty(arg)
            @warn "Ignoring malformed argument (missing '='): $arg"
        end
    end

    # Merge with defaults (user-supplied values override)
    params = merge(defaults, input_dict)

    # Convert types
    ncfile     = params["ncfile"]
    grid_size  = parse(Float64, params["grid_size"])
    crs_proj   = params["crs_proj"]
    chunk_size = parse(Int, params["chunk_size"])

    # Derived paths
    csv_path  = replace(ncfile, ".nc" => ".csv")
    meta_path = replace(ncfile, ".nc" => ".meta.json")

    # Compute CSV + metadata if missing
    if !(isfile(csv_path) && isfile(meta_path))
        println("Computing local exposure and water age...")
        compute_local_data(ncfile; grid_size=grid_size, target_crs=crs_proj, chunk_size=chunk_size)
    end

    # Export each numeric column as a single-band GeoTIFF
    println("Building geospatial files...")
    raster_paths = export_geospatial(csv_path, meta_path; fmt="GTiff")

    println("All done.")
end

# Call
main()







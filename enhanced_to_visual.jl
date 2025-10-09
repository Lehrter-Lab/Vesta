using DataFrames, GeoDataFrames, Statistics, ArchGDAL, GeometryBasics, Proj, CairoMakie, Colors, NCDatasets, ThreadsX, Base.Threads
function compute_local(ncfile::String;
                       grid_size::Float64=5_000.0,       # bin size in meters
                       target_crs::String="EPSG:5070",
                       chunk_size::Int=1_000_000)

    ds = NCDataset(ncfile, "r")
    N  = length(ds["pid"])

    ## Pass 1: determine domain bounds in target CRS
    # Temp file for bounding box checkpoint
    bbox_tmpfile = replace(ncfile, ".nc" => ".bbox.tmp")

    if isfile(bbox_tmpfile)
        println("Found bounding box checkpoint: $bbox_tmpfile")
        bbox_data = readlines(bbox_tmpfile)
        min_lon, max_lon = parse.(Float64, split(bbox_data[1]))
        min_lat, max_lat = parse.(Float64, split(bbox_data[2]))
        println("Loaded bounding box from checkpoint:")
        println("  lon: [$min_lon, $max_lon], lat: [$min_lat, $max_lat]")
    else
        println("Computing domain bounding box from file...")
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

        # Write checkpoint file (two lines: lon range and lat range)
        open(bbox_tmpfile, "w") do io
            println(io, "$min_lon $max_lon")
            println(io, "$min_lat $max_lat")
			flush(io)
        end
        println("Saved bounding box checkpoint to $bbox_tmpfile")
    end
	
	# Projection defs
	trans_fwd = Proj.Transformation("EPSG:4326", target_crs; always_xy=true)
	trans_inv = Proj.Transformation(target_crs, "EPSG:4326"; always_xy=true)
	
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

    # Build equal-area grid edges (meters in target_crs)
    edges_x = min_x:grid_size:max_x
    edges_y = min_y:grid_size:max_y
    n_x, n_y = length(edges_x)-1, length(edges_y)-1
	
    # Global accumulators + per-thread buffers
    dt_sum_cell        = zeros(Float64, n_x, n_y)
    time_weighted_cell = zeros(Float64, n_x, n_y)
    n_particles_cell   = zeros(Int, n_x, n_y)

    nthreads = Threads.nthreads()
    thread_dt = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    thread_tw = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    thread_np = [zeros(Int,     n_x, n_y) for _ in 1:nthreads]

    ## Pass 2: chunked processing
	chunk_indices = collect(1:chunk_size:N)
	for i in 1:length(chunk_indices)
		start = chunk_indices[i]
		stop  = min(start+chunk_size-1, N)

		# I/O
		@views begin
			pid_chunk  = ds["pid"][start:stop]
			lon_chunk  = ds["lon"][start:stop]
			lat_chunk  = ds["lat"][start:stop]
			time_chunk = ds["time"][start:stop]
		end

		# Transform coords to projected CRS (EPSG:5070)
        proj = trans_fwd.(lon_chunk, lat_chunk)
		if isa(proj, Tuple) && length(proj) == 2 &&
		   isa(proj[1], AbstractArray) && isa(proj[2], AbstractArray)
		    x_chunk, y_chunk = proj
		else
		    x_chunk = first.(proj)
		    y_chunk = last.(proj)
		end

		# Drop old lon/lat
        lon_chunk = nothing
        lat_chunk = nothing

		# Sort by pid
		sorted_idx  = sortperm(pid_chunk)
		pid_sorted  = pid_chunk[sorted_idx]
		x_sorted    = x_chunk[sorted_idx]
		y_sorted    = y_chunk[sorted_idx]
		time_sorted = time_chunk[sorted_idx]

		# Drop old chunk arrays
        pid_chunk  = nothing
        x_chunk    = nothing
        y_chunk    = nothing
        time_chunk = nothing

		# Find group boundaries (per-particle trajectory)
		breaks = vcat(1,
					  findall(diff(pid_sorted) .!= 0) .+ 1,
					  length(pid_sorted)+1)

		# Parallel processing over groups
		@threads for g in 1:length(breaks)-1 schedule=dynamic
			tid = threadid()
			local_dt = thread_dt[tid]
			local_tw = thread_tw[tid]
			local_np = thread_np[tid]

			lo = breaks[g]
			hi = breaks[g+1]-1

			xs = view(x_sorted, lo:hi)
			ys = view(y_sorted, lo:hi)
			ts = view(time_sorted, lo:hi)

			if length(xs) > 1
				x0, y0 = edges_x[1], edges_y[1]
				# Loop manually instead of using diff() to avoid temporary allocations
				for j in 2:length(xs)
				    dt_val = ts[j] - ts[j-1]
				    if dt_val > 0
				        x_bin = clamp(Int(floor((xs[j] - x0) / grid_size)) + 1, 1, n_x)
						y_bin = clamp(Int(floor((ys[j] - y0) / grid_size)) + 1, 1, n_y)
				        @inbounds begin
				            local_dt[x_bin, y_bin] += dt_val
				            local_tw[x_bin, y_bin] += dt_val * ts[j]
				            local_np[x_bin, y_bin] += 1
				        end
				    end
				end
			end
		end
		
		# Drop sorted arrays after processing chunk
        pid_sorted  = nothing
        x_sorted    = nothing
        y_sorted    = nothing
        time_sorted = nothing
	end

    close(ds)

    # Reduce thread-local buffers
    for tid in 1:nthreads
        dt_sum_cell        .+= thread_dt[tid]
        time_weighted_cell .+= thread_tw[tid]
        n_particles_cell   .+= thread_np[tid]
    end

	# Build DataFrame only for non-empty cells
	n_rows = count(!iszero, n_particles_cell)
	rows = Vector{NamedTuple}(undef, n_rows)

	k = 1
	for xi in 1:n_x, yi in 1:n_y
		if n_particles_cell[xi, yi] > 0
			rows[k] = (
				x_bin = xi, y_bin = yi,
				dt_sum = dt_sum_cell[xi, yi],
				time_weighted = time_weighted_cell[xi, yi],
				n_particles = n_particles_cell[xi, yi])
			k += 1
		end
	end

	df = DataFrame(rows)
    df.mean_exp_time  = df.dt_sum ./ df.n_particles
    df.mean_water_age = df.time_weighted ./ df.dt_sum

    # Create Rect geometries in target CRS (meters)
    df.geometry = [ Rect(
                        edges_x[r.x_bin],
                        edges_y[r.y_bin],
                        edges_x[r.x_bin+1] - edges_x[r.x_bin],
                        edges_y[r.y_bin+1] - edges_y[r.y_bin]
                    ) for r in eachrow(df) ]

    gdf = GeoDataFrames.GeoDataFrame(df, :geometry)
    GeoDataFrames.setcrs!(gdf, target_crs)   # already in target_crs

    output_path = replace(ncfile, ".nc" => ".gpkg")
    try
		GeoDataFrames.write(gdf, output_path)
		println("Saved local exposure (equal-area) to GeoPackage: $output_path")

		# Cleanup only after success
		if isfile(bbox_tmpfile)
			rm(bbox_tmpfile; force=true)
			println("Deleted bounding box checkpoint: $bbox_tmpfile")
		end
	catch err
		@warn "Error saving GeoPackage; bounding box checkpoint retained for resume" exception=(err, catch_backtrace())
	end
	return gdf
end
function plot_heatmap(gdf, value_col::Symbol, title::String, output_path::String; crs="EPSG:3857", cmap=:viridis, log_scale=true, units=3600)
	data = deepcopy(gdf)
	
    # convert time-like values
    if occursin("time", String(value_col)) || occursin("age", String(value_col)) || occursin("exp", String(value_col))
        data[!, value_col] ./= units
    end

    # filter log scale
    if log_scale
        data = filter(row -> row[value_col] > 0, data)
        vmin, vmax = extrema(data[!, value_col])
        norm_fn = log10
    else
        vmin, vmax = extrema(data[!, value_col])
        norm_fn = identity
    end

    # reproject geometries
    if crs != data.crs
        data = GeoDataFrames.reproject(data, crs)
    end

    polys = data.geometry
    vals  = norm_fn.(data[!, value_col])

    fig = Figure(resolution=(1200,800))
    ax  = Axis(fig[1,1], title=title)
    poly!(ax, polys, color=vals, colormap=cmap)
    Colorbar(fig[1,2], limits=(vmin,vmax), colormap=cmap, label=title)
    fig[0,:] = Label(fig, title, fontsize=18)

    save(output_path, fig)
    return fig
end
function main()
    # Input parameters ------------------------------------------------
    ncfile       = "particle_enhanced.nc"
    grid_size    = 5000.0       # meters
    crs_proj     = "EPSG:5070"  # equal-area projection
    chunk_size   = 1000000
    sample_N     = 10000        # particles per Monte Carlo iteration
    iterations_M = 50           # Monte Carlo iterations

    # Compute local exposure and water age --------------------------------
    println("Computing local exposure and water age...")
    gdf = compute_local(ncfile;
                        grid_size=grid_size,
                        target_crs=crs_proj,
                        chunk_size=chunk_size)

    println("Done computing local values.")

    # Plot heatmap of mean exposure time ---------------------------------
    println("Plotting heatmap...")
    plot_heatmap(gdf, :mean_exp_time, "Mean Exposure Time", "heatmap_exposure.png"; crs="EPSG:3857")

    println("All done.")
end
main()



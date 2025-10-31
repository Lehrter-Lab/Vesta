using DataFrames, CSV, Statistics, JSON3
using ArchGDAL, GeometryBasics, GeoDataFrames, Proj, NCDatasets, GeoTables
using CairoMakie, Colors, ThreadsX, Base.Threads
using FilePathsBase

# Compute and Save Tabular Results (CSV + Metadata)
function compute_local_data(ncfile::String;
                            grid_size::Float64=5_000.0,
                            target_crs::String="EPSG:5070",
                            chunk_size::Int=1_000_000)

    ds = NCDataset(ncfile, "r")
    N  = length(ds["pid"])

    # Bounding Box Checkpoint Logic
    bbox_tmpfile = replace(ncfile, ".nc" => ".bbox.tmp")

    if isfile(bbox_tmpfile)
        println("Found bounding box checkpoint: $bbox_tmpfile")
        bbox_data = readlines(bbox_tmpfile)
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

    # Projection Setup
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
    edges_x = min_x:grid_size:max_x
    edges_y = min_y:grid_size:max_y
    n_x, n_y = length(edges_x)-1, length(edges_y)-1

    # Initialize Buffers
    dt_sum_cell        = zeros(Float64, n_x, n_y)
    time_weighted_cell = zeros(Float64, n_x, n_y)
    n_particles_cell   = zeros(Int, n_x, n_y)

    nthreads = Threads.nthreads()
    thread_dt = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    thread_tw = [zeros(Float64, n_x, n_y) for _ in 1:nthreads]
    thread_np = [zeros(Int,     n_x, n_y) for _ in 1:nthreads]

    println("Pass 2: Processing in chunks...")

    chunk_indices = collect(1:chunk_size:N)
    for i in 1:length(chunk_indices)
        start = chunk_indices[i]
        stop  = min(start+chunk_size-1, N)

        @views begin
            pid_chunk  = ds["pid"][start:stop]
            lon_chunk  = ds["lon"][start:stop]
            lat_chunk  = ds["lat"][start:stop]
            time_chunk = ds["time"][start:stop]
        end

        # Transform coordinates (robust version)
        proj_result = trans_fwd.(lon_chunk, lat_chunk)
        if isa(proj_result, Tuple) && length(proj_result) == 2 &&
           isa(proj_result[1], AbstractArray) && isa(proj_result[2], AbstractArray)
            x_chunk, y_chunk = proj_result
        else
            x_chunk = first.(proj_result)
            y_chunk = last.(proj_result)
        end

        lon_chunk = nothing
        lat_chunk = nothing

        # Sort by particle id
        sorted_idx  = sortperm(pid_chunk)
        pid_sorted  = pid_chunk[sorted_idx]
        x_sorted    = x_chunk[sorted_idx]
        y_sorted    = y_chunk[sorted_idx]
        time_sorted = time_chunk[sorted_idx]

        pid_chunk  = nothing
        x_chunk    = nothing
        y_chunk    = nothing
        time_chunk = nothing

        # Group by pid
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
    end

    close(ds)

    for tid in 1:nthreads
        dt_sum_cell        .+= thread_dt[tid]
        time_weighted_cell .+= thread_tw[tid]
        n_particles_cell   .+= thread_np[tid]
    end

    n_rows = count(!iszero, n_particles_cell)
    rows = Vector{NamedTuple}(undef, n_rows)
    k = 1
    for xi in 1:n_x, yi in 1:n_y
        if n_particles_cell[xi, yi] > 0
            rows[k] = (
                x_bin = xi, y_bin = yi,
                dt_sum = dt_sum_cell[xi, yi],
                time_weighted = time_weighted_cell[xi, yi],
                n_particles = n_particles_cell[xi, yi]
            )
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
    meta = Dict(
        "grid_size" => grid_size,
        "target_crs" => target_crs,
        "bounds" => (min_x, max_x, min_y, max_y),
        "n_x" => n_x,
        "n_y" => n_y
    )
    JSON3.write(meta_path, meta)
    println("Saved compact metadata to $meta_path")

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

function export_geospatial(csv_path::String, meta_path::String; fmt::String="GPKG")
    using CSV, DataFrames, GeometryBasics, GeoTables, ArchGDAL, JSON3

    df = CSV.read(csv_path, DataFrame)
    meta = JSON3.read(Base.read(meta_path, String))

    grid_size = meta["grid_size"]
    (min_x, max_x, min_y, max_y) = meta["bounds"]
    n_x, n_y = meta["n_x"], meta["n_y"]
    target_crs = meta["target_crs"]

    edges_x = collect(min_x:grid_size:max_x)
    edges_y = collect(min_y:grid_size:max_y)

    df.geometry = [
        GeometryBasics.Polygon(Point2f[
            Point2f(edges_x[r.x_bin],     edges_y[r.y_bin]),
            Point2f(edges_x[r.x_bin+1],   edges_y[r.y_bin]),
            Point2f(edges_x[r.x_bin+1],   edges_y[r.y_bin+1]),
            Point2f(edges_x[r.x_bin],     edges_y[r.y_bin+1]),
            Point2f(edges_x[r.x_bin],     edges_y[r.y_bin])
        ])
        for r in eachrow(df)
    ]

    # Replace GeoDataFrame with GeoTable
    gdf = GeoTables.GeoTable(df; geometry=:geometry, crs=target_crs)
    output_path = replace(csv_path, ".csv" => fmt == "SHP" ? ".shp" : ".gpkg")

    ArchGDAL.create(output_path, driver=fmt) do dataset
        srs = ArchGDAL.importEPSG(parse(Int, split(target_crs, ":")[2]))
        layer = ArchGDAL.create_layer(dataset, "grid_data"; srs=srs)

        for c in names(df)
            if c != :geometry
                t = eltype(df[!, c]) <: Integer ? ArchGDAL.OFTInteger : ArchGDAL.OFTReal
                ArchGDAL.createfield(layer, String(c) => t)
            end
        end

        for r in eachrow(df)
            geom = ArchGDAL.creategeometry(r.geometry)
            attrs = Dict(Symbol(k)=>r[k] for k in names(df) if k != :geometry)
            ArchGDAL.createfeature(layer, geom, attrs)
        end
    end

    println("Exported geospatial data to $output_path")
    return gdf
end

# Plot Heatmap
function plot_heatmap(gdf, value_col::Symbol, title::String, output_path::String;
                      crs="EPSG:3857", cmap=:viridis, log_scale=true, units=3600)
    data = deepcopy(gdf)

    if occursin("time", String(value_col)) || occursin("age", String(value_col)) || occursin("exp", String(value_col))
        data[!, value_col] ./= units
    end

    if log_scale
        data = filter(row -> row[value_col] > 0, data)
        if isempty(data)  # change #7
            @warn "No valid data for plotting $value_col"
            return nothing
        end
        vmin, vmax = extrema(data[!, value_col])
        norm_fn = log10
    else
        vmin, vmax = extrema(data[!, value_col])
        norm_fn = identity
    end

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
    println("Saved plot to $output_path")
    return fig
end

# Main
function main(; resume=false)
    ncfile = "particle_enhanced.nc"
    grid_size = 5000.0
    crs_proj = "EPSG:5070"
    chunk_size = 10_000_000

    csv_path  = replace(ncfile, ".nc" => ".csv")
    meta_path = replace(ncfile, ".nc" => ".meta.json")

    # Check for existing CSV + metadata
    if isfile(csv_path) && isfile(meta_path)
        println("Detected existing CSV and metadata files.")
        println("Skipping computation and proceeding to geospatial build...")
    else
        println("Computing local exposure and water age...")
        compute_local_data(ncfile; grid_size=grid_size, target_crs=crs_proj, chunk_size=chunk_size)
    end

    # Proceed to geospatial export 
    println("Building geospatial file...")
    gdf = export_geospatial(csv_path, meta_path; fmt="GPKG")

    # Plot results
    println("Plotting heatmap...")
    plot_heatmap(gdf, :mean_exp_time, "Mean Exposure Time", "heatmap_exposure.png"; crs="EPSG:3857")

    println("All done.")
end

# Call
main()

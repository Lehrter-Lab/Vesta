using ZipFile, Shapefile, ArchGDAL, GeometryBasics, NCDatasets, Dates, Printf, Base.Threads
function unzip_shapefile(zip_path::String, extract_dir::String)
    isdir(extract_dir) || mkdir(extract_dir)
    archive = ZipFile.Reader(zip_path)
    for file in archive.files
        open(joinpath(extract_dir, file.name), "w") do io
            write(io, read(file))
        end
    end
    close(archive)
    for f in readdir(extract_dir)
        endswith(f, ".shp") && return joinpath(extract_dir, f)
    end
    error("No .shp file found in archive")
end
function load_domain_geometry(shp_path::String)
    shp = Shapefile.Table(shp_path)
    arch_polygons = Vector{ArchGDAL.IGeometry}(undef, length(Shapefile.shapes(shp)))
    idx = 1

    for poly in Shapefile.shapes(shp)
        coords = hasproperty(poly, :coordinates) ? coordinates(poly) :
                 hasproperty(poly, :points)      ? [(pt.x, pt.y) for pt in poly.points] :
                 error("Cannot extract coordinates from polygon of type $(typeof(poly))")

        rings = isa(coords[1], Tuple) ? [coords] : coords
        for ring in rings
            ring = ring[1] == ring[end] ? ring : vcat(ring, ring[1])
            arch_polygons[idx] = ArchGDAL.createpolygon(ring)
            idx += 1
        end
    end

    return reduce((a,b) -> ArchGDAL.union(a,b), arch_polygons[1:idx-1])
end
function compute_inside_flags(lons::Vector{Float32}, lats::Vector{Float32}, domain_geom::ArchGDAL.IGeometry)
    n = length(lons)
    inside_flags = falses(n)

    # Bounding box pre-filter
    env = ArchGDAL.envelope(domain_geom)
    minx, maxx = env.MinX, env.MaxX
    miny, maxy = env.MinY, env.MaxY
    idxs = findall((lons .>= minx) .& (lons .<= maxx) .& (lats .>= miny) .& (lats .<= maxy))

    # Threaded inside-domain check
    Threads.@threads for i in idxs
        pt = ArchGDAL.createpoint(lons[i], lats[i])
        inside_flags[i] = ArchGDAL.contains(domain_geom, pt)
    end

    return inside_flags
end
function save_particle_enhanced_parallel(path_nc::String, output_nc::String, domain_geom::ArchGDAL.IGeometry; chunk_size::Int=1_000_000)
    ds = NCDataset(path_nc, "r")
    N = length(ds["pid"][:])
    inside_vec = falses(N)

    n_chunks = ceil(Int, N / chunk_size)

    Threads.@threads for c in 1:n_chunks
        start = (c - 1) * chunk_size + 1
        stop = min(c * chunk_size, N)
        lon_chunk = ds["lon"][start:stop]
        lat_chunk = ds["lat"][start:stop]

        # Compute inside flags for this chunk
        inside_chunk = compute_inside_flags(lon_chunk, lat_chunk, domain_geom)
        inside_vec[start:stop] .= inside_chunk  # in-place assignment
    end

    # Save enhanced NetCDF
    ds_out = NCDataset(output_nc, "c")
    defDim(ds_out, "point", N)
    var_pid = defVar(ds_out, "pid", Int32, ("point",))
    var_lon = defVar(ds_out, "lon", Float32, ("point",))
    var_lat = defVar(ds_out, "lat", Float32, ("point",))
    var_time = defVar(ds_out, "time", Float32, ("point",))
    var_inside = defVar(ds_out, "inside", UInt8, ("point",))

    var_pid[:] = ds["pid"][:]
    var_lon[:] = ds["lon"][:]
    var_lat[:] = ds["lat"][:]
    var_time[:] = ds["time"][:]
    var_inside[:] = UInt8.(inside_vec)

    close(ds)
    close(ds_out)
end
function compute_times_chunked(input_nc::String, output_nc::String; chunk_size::Int=1_000_000)
    ds = NCDataset(input_nc, "r")
    Npoints = length(ds["pid"][:])
    all_pids = ds["pid"][:]
    max_pid = maximum(all_pids)

    # Output file
    ds_out = NCDataset(output_nc, "c")
    defDim(ds_out, "pid", max_pid)
    var_res_time = defVar(ds_out, "res_time", Float32, ("pid",))
    var_exp_time = defVar(ds_out, "exp_time", Float32, ("pid",))

    # Global state, tracked per particle ID
    first_exit_found = falses(max_pid)
    res_times   = zeros(Float32, max_pid)
    exp_times   = zeros(Float32, max_pid)
    last_time   = fill(-1.0, max_pid)  # sentinel meaning "no time seen yet"

    n_chunks = ceil(Int, Npoints / chunk_size)
    println("Processing $n_chunks chunks for $Npoints points across $(max_pid) pids...")

    Threads.@threads for c in 1:n_chunks
        start_idx = (c - 1) * chunk_size + 1
        stop_idx  = min(c * chunk_size, Npoints)
        chunk_idx = start_idx:stop_idx

        inside = ds["inside"][chunk_idx]
        times  = ds["time"][chunk_idx]
        pids   = all_pids[chunk_idx]

        @inbounds for i_local in eachindex(chunk_idx)
            pid = pids[i_local]
            t   = times[i_local]

            # Residence time
            if !first_exit_found[pid]
                if inside[i_local] == 0  # exited here
                    if last_time[pid] > 0.0
                        res_times[pid] = t - last_time[pid]  # residence = exit - entry
                        first_exit_found[pid] = true
                    end
                end
            end

            # Exposure time (accumulate delta_T while inside)
            if last_time[pid] > 0.0
                delta_T = t - last_time[pid]
                if inside[i_local] == 1
                    exp_times[pid] += delta_T
                end
            end

            # Update last_time always
            last_time[pid] = t
        end
    end

    #Write results
    var_res_time[:] = res_times
    var_exp_time[:] = exp_times

    close(ds)
    close(ds_out)
    println("Residence and exposure times saved to $output_nc")
end
function main()
    zip_shp     = "../../juliaParticle/bbox_dissolve.zip"      # your shapefile zip
    extract_dir = "./shapefile_extracted"  # temp folder for unzip
    input_nc    = "./particleFall.nc"           # input particle NetCDF
    enhanced_nc = "./particle_enhanced.nc"  # output NetCDF with 'inside' flag
    times_nc    = "./particle_times.nc"     # output NetCDF with res/exp times

    shp_path = unzip_shapefile(zip_shp, extract_dir)
    println("Shapefile extracted to: $shp_path")
    domain_geom = load_domain_geometry(shp_path)
    println("Domain geometry loaded.")

    println("Computing inside flags...")
    save_particle_enhanced_parallel(input_nc, enhanced_nc, domain_geom; chunk_size=1_000_000)
    println("Enhanced NetCDF saved to $enhanced_nc")

    println("Computing residence/exposure times...")
    compute_times_chunked(enhanced_nc, times_nc; chunk_size=1_000_000)
    println("Times NetCDF saved to $times_nc")
end
main()

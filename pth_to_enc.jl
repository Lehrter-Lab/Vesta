using Shapefile, ZipFile, ArchGDAL, Proj4, NCDatasets, Printf, Base.Threads

# Unzip shapefile and get .shp + .prj
function unzip_shapefile(zip_path::String, extract_dir::String)
    isdir(extract_dir) || mkpath(extract_dir)
    shp_path, prj_text = nothing, nothing
    archive = ZipFile.Reader(zip_path)
    for f in archive.files
        outpath = joinpath(extract_dir, f.name)
        mkpath(dirname(outpath))
        open(outpath, "w") do io
            write(io, read(f))
        end
        if endswith(lowercase(f.name), ".shp"); shp_path = outpath end
        if endswith(lowercase(f.name), ".prj"); prj_text = read(outpath, String) end
    end
    close(archive)
    shp_path === nothing && error("No .shp found in $zip_path")
    return shp_path, prj_text
end

# Load polygons from shapefile
function load_polygons(shp_path::String)
    table = Shapefile.Table(shp_path)
    polygons = Vector{Vector{Vector{Tuple{Float64,Float64}}}}()
    for shape in Shapefile.shapes(table)
        coords = hasproperty(shape, :coordinates) ? coordinates(shape) :
                 hasproperty(shape, :points)      ? [(pt.x, pt.y) for pt in shape.points] :
                 error("Unknown polygon type")
        rings = isa(coords[1], Tuple) ? [coords] : coords
        push!(polygons, [ [ (Float64(p[1]), Float64(p[2])) for p in ring ] for ring in rings ])
    end
    return polygons
end

# Parse EPSG from .prj text
function parse_epsg(prj_text::Union{String,Nothing})
    if prj_text === nothing; return nothing end
    m = match(Regex("EPSG:?[0-9]{3,6}", "i"), prj_text)  # "i" for ignore case
    if m !== nothing
        # Extract the digits using a capturing group
        m_digits = match(r"[0-9]{3,6}", m.match)
        return parse(Int, m_digits.match)
    else
        return nothing
    end
    return m !== nothing ? parse(Int, m.captures[1]) : nothing
end

# Reproject polygons to NAD83 (EPSG:4269)
function reproject_polygons(polygons::Vector{Vector{Vector{Tuple{Float64,Float64}}}}, src_epsg::Int)
    # Set up transformation from source EPSG to NAD83
    src_proj = Proj4.Transformation("EPSG:$src_epsg", "EPSG:4269")

    # Loop over each polygon and ring, transforming each point
    reprojected = Vector{Vector{Vector{Tuple{Float64,Float64}}}}(undef, length(polygons))
    for i_poly in 1:length(polygons)
        poly = polygons[i_poly]
        reprojected[i_poly] = Vector{Vector{Tuple{Float64,Float64}}}(undef, length(poly))
        for i_ring in 1:length(poly)
            ring = poly[i_ring]
            reprojected[i_poly][i_ring] = [(Proj4.transform(src_proj, p[1], p[2])[1],
                                            Proj4.transform(src_proj, p[1], p[2])[2]) for p in ring]
        end
    end
    return reprojected
end

# Combine polygons into single ArchGDAL geometry
function make_domain_geom(polygons)
    geoms = ArchGDAL.IGeometry[]
    for poly in polygons, ring in poly
        ring_closed = ring[1] == ring[end] ? ring : vcat(ring, ring[1])
        push!(geoms, ArchGDAL.createpolygon(ring_closed))
    end
    g = geoms[1]
    for i in 2:length(geoms); g = ArchGDAL.union(g, geoms[i]) end
    return g
end

# Compute inside flags for particles
function compute_inside(lons::Vector{Float64}, lats::Vector{Float64}, domain_geom::ArchGDAL.IGeometry)
    n      = length(lons)
    inside = falses(n)
    env    = ArchGDAL.envelope(domain_geom)
    idxs   = findall((lons .>= env.MinX) .& (lons .<= env.MaxX) .& (lats .>= env.MinY) .& (lats .<= env.MaxY))
    for i in idxs
        pt = ArchGDAL.createpoint(lons[i], lats[i])
        inside[i] = ArchGDAL.contains(domain_geom, pt)
    end
    return inside
end

# Initialize enhanced NetCDF
function init_enhanced_nc(file::String; compress=4)
    ds = NCDataset(file, "c")
    defDim(ds, "point", Inf)
    v_pid    = defVar(ds, "pid", Int32, ("point",))
    v_lon    = defVar(ds, "lon", Float32, ("point",); deflatelevel=compress, shuffle=true)
    v_lat    = defVar(ds, "lat", Float32, ("point",); deflatelevel=compress, shuffle=true)
    v_depth  = defVar(ds, "depth", Float32, ("point",); deflatelevel=compress, shuffle=true)
    v_time   = defVar(ds, "time", Float64, ("point",); deflatelevel=compress, shuffle=true)
    v_inside = defVar(ds, "inside", UInt8, ("point",); deflatelevel=compress, shuffle=true)
    return ds, (v_pid,v_lon,v_lat,v_depth,v_time,v_inside)
end
# Append a timestep to NetCDF
function append_timestep!(ds::NCDataset, vars, pidv, lonv, latv, depthv, timev, insidev, start_idx::Int)
    n = length(pidv)
    end_idx = start_idx + n - 1
    v_pid,v_lon,v_lat,v_depth,v_time,v_inside = vars
    v_pid[start_idx:end_idx]    = Int32.(pidv)
    v_lon[start_idx:end_idx]    = Float32.(lonv)
    v_lat[start_idx:end_idx]    = Float32.(latv)
    v_depth[start_idx:end_idx]  = Float32.(depthv)
    v_time[start_idx:end_idx]   = Float64.(timev)
    v_inside[start_idx:end_idx] = UInt8.(insidev)
    return end_idx + 1
end

# Compute residence/exposure times
function compute_times(enhanced_nc::String, output_nc::String)
    ds         = NCDataset(enhanced_nc,"r")
    pids       = Int.(ds["pid"][:])
    inside     = ds["inside"][:]
    times      = ds["time"][:]
    max_pid    = maximum(pids)
    res_time   = zeros(Float32, max_pid)
    exp_time   = zeros(Float32, max_pid)
    entry      = fill(NaN,max_pid)
    last       = fill(NaN,max_pid)
    first_exit = falses(max_pid)
    ds_out     = NCDataset(output_nc,"c")
    defDim(ds_out,"pid",max_pid)
    var_res = defVar(ds_out,"res_time",Float32,("pid",))
    var_exp = defVar(ds_out,"exp_time",Float32,("pid",))
    for i in eachindex(pids)
        pid = pids[i]; t = times[i]; inside_flag = inside[i]==1
        if isnan(entry[pid]) && inside_flag; entry[pid]=t end
        # Residence time
        if !first_exit[pid] && !inside_flag && !isnan(entry[pid])
            res_time[pid]   = Float32(t-entry[pid])
            first_exit[pid] = true
        end
        # Exposure time
        if !isnan(last[pid])
            dt = t-last[pid]
            if inside_flag; exp_time[pid]+=Float32(dt) end
        end
        last[pid]=t
    end
    var_res[:] = res_time
    var_exp[:] = exp_time
    close(ds); close(ds_out)
end

# Main streaming pipeline
function pth_to_enhanced_main(pth_file::String, shp_zip::String,
                              enhanced_out::String, times_out::String)
    shp_path, prj_text = unzip_shapefile(shp_zip, "./_shp_extracted")
    polygons = load_polygons(shp_path)
    src_epsg = parse_epsg(prj_text)
    if src_epsg !== nothing && src_epsg != 4269
        @warn "Shapefile EPSG $src_epsg -> reprojecting to NAD83 (4269)"
        polygons = reproject_polygons(polygons, src_epsg)
    elseif src_epsg===nothing
        @warn "Shapefile EPSG not detected; assuming NAD83 (4269)"
    end
    domain_geom  = make_domain_geom(polygons)
    ds_enh, vars = init_enhanced_nc(enhanced_out)
    next_idx = 1; max_pid = 0
    open(pth_file,"r") do io
        while !eof(io)
            header = readline(io)
            parts  = split(strip(header))
            length(parts)!=2 && continue
            time_sec, n_particles = try parse(Float64, parts[1]), parse(Int, parts[2]) catch; continue end
            pidv   = Vector{Int}(undef,n_particles)
            lonv   = Vector{Float64}(undef,n_particles)
            latv   = Vector{Float64}(undef,n_particles)
            depthv = Vector{Float64}(undef,n_particles)
            timev  = Vector{Float64}(undef,n_particles)
            for i in 1:n_particles
                line=readline(io)
                vals=split(strip(line))
                length(vals)<4 && error("Malformed line: $line")
                pidv[i]=parse(Int,vals[1])
                lonv[i]=parse(Float64,vals[2])
                latv[i]=parse(Float64,vals[3])
                depthv[i]=parse(Float64,vals[4])
                timev[i]=time_sec
                max_pid = max(max_pid,pidv[i])
            end
            inside_flags = compute_inside(lonv, latv, domain_geom)
            next_idx = append_timestep!(ds_enh, vars, pidv, lonv, latv, depthv, timev, inside_flags, next_idx)
            @info "Wrote timestep $(time_sec) n_particles=$(n_particles)"
        end
    end
    close(ds_enh)
    @info "Enhanced NetCDF done: $enhanced_out"
    compute_times(enhanced_out, times_out)
    @info "Residence/exposure times done: $times_out"
end

# CLI interface with defaults
function main(; resume=false)
    # Example STDIN input:
    # echo "pth_file=particle.pth shp_zip=../../juliaParticle/bbox_dissolve.zip enhanced_out=particle_enhanced.nc times_out=particle_times.nc" | julia pth_to_enhanced.jl

    # Allowed keys and defaults
    valid_keys = Set(["pth_file", "shp_zip", "enhanced_out", "times_out"])
    defaults = Dict("pth_file"     => "particle.pth",
                    "shp_zip"      => "../../juliaParticle/bbox_dissolve.zip",
                    "enhanced_out" => "particle_enhanced.nc",
                    "times_out"    => "particle_times.nc")

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

    # Call main processing function
    pth_file     = params["pth_file"]
    shp_zip      = params["shp_zip"]
    enhanced_out = params["enhanced_out"]
    times_out    = params["times_out"]

    pth_to_enhanced_main(pth_file, shp_zip, enhanced_out, times_out)
end

# Call
main()

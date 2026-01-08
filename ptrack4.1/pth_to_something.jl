# Use hgrid subset to find bool for inside relevant domain
function compute_inside(elements,ingrid)
	lines = readlines(ingrid)
    # Second line has ne np
    ne, np = parse.(Int, split(lines[2]))
    # Element lines start after the nodes
    ele_lines = lines[2 + np + 1 : 2 + np + ne]
    # Extract element IDs
    ele = [parse(Int, split(line)[1]) for line in ele_lines]
	subdomain = Set(ele)
    return in.(elements, Ref(subdomain))
	
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
function append_timestep!(ds::NCDataset, vars, pidv, elev, lonv, latv, depthv, timev, insidev, start_idx::Int)
    n = length(pidv)
    end_idx = start_idx + n - 1
    v_pid,v_lon,v_lat,v_depth,v_time,v_inside = vars
    v_pid[start_idx:end_idx]    = Int32.(pidv)
	v_ele[start_idx:end_idx]    = Int32.(elev)	
    v_lon[start_idx:end_idx]    = Float32.(lonv)
    v_lat[start_idx:end_idx]    = Float32.(latv)
    v_depth[start_idx:end_idx]  = Float32.(depthv)
    v_time[start_idx:end_idx]   = Float64.(timev)
    v_inside[start_idx:end_idx] = UInt8.(insidev)
    return end_idx + 1
end
	
# Read .pth file
open(pth_file,"r") do io
        while !eof(io)
			# Handle header
            header = readline(io)
            parts  = split(strip(header))
            length(parts)!=2 && continue
            time_sec, n_particles = try parse(Float64, parts[1]), parse(Int, parts[2]) catch; continue end
			# Initialize empty vectors for number of particles
            pidv   = Vector{Int}(undef,n_particles)
	        elev   = Vector{Int}(undef,n_particles)
            lonv   = Vector{Float64}(undef,n_particles)
            latv   = Vector{Float64}(undef,n_particles)
            depthv = Vector{Float64}(undef,n_particles)
            timev  = Vector{Float64}(undef,n_particles)
			# Fill
            for i in 1:n_particles
                line=readline(io)
                vals=split(strip(line))
                length(vals)<5 && error("Malformed line: $line")
                pidv[i]   = parse(Int,vals[1])
		        elev[i]   = parse(Int,vals[2])
                lonv[i]   = parse(Float64,vals[3])
                latv[i]   = parse(Float64,vals[4])
                depthv[i] = parse(Float64,vals[5])
                timev[i]  = time_sec
                max_pid   = max(max_pid,pidv[i])
            end
            inside_flags = compute_inside(elev, domain_geom)
            next_idx     = append_timestep!(ds_enh, vars, pidv, elev, lonv, latv, depthv, timev, inside_flags, next_idx)
            @info "Wrote timestep $(time_sec) n_particles=$(n_particles)"
        end
    end
    close(ds_enh)
    @info "Enhanced NetCDF done: $enhanced_out"
end

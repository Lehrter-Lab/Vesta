function subdomain_gr3(infile)
	lines = readlines(filename)
    # Second line has ne np
    ne, np = parse.(Int, split(lines[2]))
    # Element lines start after the nodes
    ele_lines = lines[2 + np + 1 : 2 + np + ne]
    # Extract element IDs
    ele = [parse(Int, split(line)[1]) for line in ele_lines]
    return ele

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
            inside_flags = compute_inside(lonv, latv, domain_geom)
            next_idx     = append_timestep!(ds_enh, vars, pidv, lonv, latv, depthv, timev, inside_flags, next_idx)
            @info "Wrote timestep $(time_sec) n_particles=$(n_particles)"
        end
    end
    close(ds_enh)
    @info "Enhanced NetCDF done: $enhanced_out"
end

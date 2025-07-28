import geopandas as gpd
import numpy as np
import h5py
import os

## Input params ---------------------------------------------------------------
# Files
infile      = "olivine.h5"
output_file = "particles_new.bp"
# Times
start_time     = 2592000     # start in seconds per model record
time_interval  = 7200        # time between injections
n_time_steps   = 24          # number of injections
times          = [start_time + i * time_interval for i in range(n_time_steps)]
# .bp params
h0       = 0.01       # initial water depth (m)
rnday    = 90.0       # total model run duration (days)
dtm      = 50.0       # model timestep (seconds)
nspool   = 72         # output frequency: number of model steps between outputs
ihfskip  = 1728       # number of output steps to skip at start (spin-up)
ndeltp   = 40         # internal particle time steps per model timestep

## Param sanity check ---------------------------------------------------------
seconds_per_hour = 3600
if abs(dtm * nspool - seconds_per_hour) > 1e-6:
    raise ValueError(f"Sanity check failed: dtm * nspool = {dtm * nspool} "
                     f"but should equal {seconds_per_hour} (one hour in seconds).")

## Core funcs -----------------------------------------------------------------
# Read points from either shapefile or HDF5
def read_points(infile):
    ext = os.path.splitext(infile)[1].lower()

    if ext == ".shp":
        gdf = gpd.read_file(infile)
        if not all(gdf.geometry.geom_type == 'Point'):
            raise ValueError("Shapefile must contain only Point geometries.")
        xs = gdf.geometry.x.tolist()
        ys = gdf.geometry.y.tolist()
        if "elevation" in gdf.columns:
            zs = gdf["elevation"].tolist()
        else:
            zs = [0.0] * len(xs)
        return xs, ys, zs

    elif ext == ".h5":
        with h5py.File(infile, 'r') as f:
            node_path = None
            z_path = None

            def find_datasets(name, obj):
                nonlocal node_path, z_path
                if isinstance(obj, h5py.Dataset):
                    if name.endswith("Nodes/NodeLocs"):
                        node_path = name
                    elif name.endswith("Datasets/Z/Values"):
                        z_path = name

            f.visititems(find_datasets)

            if node_path is None or z_path is None:
                raise KeyError("Required dataset(s) not found in input file.")

            coords = f[node_path][:]  # shape (N, 3)
            xs = coords[:, 0]
            ys = coords[:, 1]
            zs = f[z_path][0, :]  # shape (1, N) → flatten to (N,)
            
            if len(zs) != len(xs):
                raise ValueError("Mismatch between number of coordinate points and z-values.")

            return xs.tolist(), ys.tolist(), zs.tolist()
    else:
        raise ValueError(f"Unsupported input file format: {ext}")
        
# Generate properly formatted particle.bp block
def generate_particle_file(xs, ys, zs, times, output_file,
                           dz=0.5,
                           nscreen=0,
                           mod_part=0,
                           ibiofoul=0,
                           ibf=1,
                           istiff=0,
                           ibnd_beh=0,
                           ics_block="2 -122.6 37.38 ics slam0 sfea0",
                           time_params="0.01 90. 50. 72 1728 40"):
    if not (len(xs) == len(ys) == len(zs)):
        raise ValueError("Length mismatch among xs, ys, zs.")
    
    particle_lines = []
    particle_id    = 1
 
    for lon, lat, max_depth in zip(xs, ys, zs):
        max_depth    = abs(float(max_depth))  # ensure positive
        local_depths = [round(-d, 2) for d in np.arange(0, max_depth + dz, dz)]

        for time in times:
            for depth in local_depths:
                line = f"{particle_id}\t{time}\t{lon:.7f}\t{lat:.7f}\t{depth}"
                particle_lines.append(line)
                particle_id += 1

    nparticles = len(particle_lines)

    # Header
    header_lines = [
        "Input for ptrack*",
        f"{nscreen} nscreen",
        f"{mod_part} mod_part (0: passive; 1: oil spill)",
        f"{ibiofoul} ibiofoul (0: no biofouling; 1: with biofouling)",
        f"{ibf} ibf !(1: forward; -1: backward)",
        f"{istiff} istiff !1: fixed distance from surface",
        f"{ibnd_beh} ibnd_beh (behavior near bnd or wet/dry; 0: reflect off; 1: slide)",
        f"{ics_block} (from param.nml)",
        f"{time_params} !h0,rnday,dtm,nspool,ihfskip,ndeltp",
        f"{nparticles} !# of particles"
    ]

    # Footer
    footer_lines = [
        "!Oil spill parameters needed only if mod_part=1",
        "1  3.0  0.2           !ihdf(0: constant diffusivity as hdc; 1: Smagorinsky),hdc,horcon",
        "1  0                  !ibuoy (buoyancy off (0) or on (1)),iwind (wind effect off (0) or on (1))",
        "20.0                  !set minimum percentage of stranding on shore (particles may be stranded if the random # exceeds this threshold)",
        "!start of biofouling paramaters, needed only if ibiofoul=1",
        "2.5e-3 0. 6.e-6 0.83 1.38   !bio_R0,bio_BT0,bio_BR,bio_den0 (ρ₀),bio_den_biolayer (ρ_D)"
    ]

    # Write to file
    with open(output_file, 'w') as f:
        f.write("\n".join(header_lines) + "\n")
        f.write("\n".join(particle_lines) + "\n")
        f.write("\n".join(footer_lines) + "\n")

    print(f"Generated {nparticles} particles and wrote to {output_file}")

## Do the work ----------------------------------------------------------------
xs, ys, zs = read_points(infile)
generate_particle_file(xs, ys, zs, times, output_file)

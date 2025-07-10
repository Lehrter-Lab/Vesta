import geopandas as gpd

# Input parameters
shapefile   = ""
output_file = "particles_new.bp"
times = [2592000, 2599200]
depths = [round(-0.5 * i, 1) for i in range(16)]

# Read a shapefile of points
def read_shapefile_points(shapefile_path):
    gdf = gpd.read_file(shapefile_path)

    if not all(gdf.geometry.geom_type == 'Point'):
        raise ValueError("Shapefile must contain only Point geometries.")

    xs = gdf.geometry.x.tolist()
    ys = gdf.geometry.y.tolist()
    return xs, ys

# Generate properly formatted particle.bp block
def generate_particle_file(xs, ys, times, depths, output_file):
    if len(xs) != len(ys):
        raise ValueError("Longitude and latitude lists must be of the same length.")

    particle_lines = []
    particle_id = 1

    for lon, lat in zip(xs, ys):
        for time in times:
            for depth in depths:
                line = f"{particle_id}\t{time}\t{lon:.7f}\t{lat:.7f}\t{depth}"
                particle_lines.append(line)
                particle_id += 1

    # Write output file
    with open(output_file, 'w') as f:
        f.write(f"{len(particle_lines)} !# of particles\n")
        f.write("\n".join(particle_lines))
    print(f"Generated {len(particle_lines)} particles and wrote to {output_file}")

## Do the work
#xs, ys = read_shapefile_points(shapefile)
xs, ys = [-70.0],[40.0]
generate_particle_file(xs, ys, times, depths, output_file)
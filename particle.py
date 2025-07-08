# CMikolaitis @ USA/DISL

import pandas as pd
import geopandas as gpd
import re

# Input files
bp     = 'particle.bp'
infile = 'particle.pth'

# Determine timestep from particle.bp
with open("particle.bp") as f:
    line9 = f.readlines()[8]
dtm, nspool = map(float, line9.strip().split()[2:4])
timestep = dtm * nspool  # seconds per record

# Load particle data
data = []
with open(infile) as f:
    for line in f:
        if m := re.match(r"\s*(\d+\.\d+)\s+\d+", line):
            t = float(m[1])
        elif m := re.match(r"\s*(\d+)\s+([\dE+-.]+)\s+([\dE+-.]+)\s+([\dE+-.]+)", line):
            data.append([t, int(m[1]), float(m[2]), float(m[3]), float(m[4])])
df = pd.DataFrame(data, columns=["time", "index", "x", "y", "z"])

# Convert particle data to geodataframe
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x'], df['y']), crs="EPSG:4326")

# Load polygon and convert to geoseries
bbox = gpd.read_file("zip:bbox_dissolve.zip").to_crs("EPSG:4326")
poly = bbox.geometry.union_all()

# Are the particles inside the polygon
gdf["inside"] = gdf.geometry.within(poly)

# Compute residence and exposure times
results = []
for idx, group in gdf.groupby("index"):
    inside = group["inside"].values
    durations = [timestep] * (len(inside) - 1) + [0]

    res_time = 0
    exp_time = 0
    current_res = 0

    for i in range(len(inside)):
        if inside[i]:
            exp_time += durations[i]
            current_res += durations[i]
        else:
            res_time = max(res_time, current_res)
            current_res = 0
    res_time = max(res_time, current_res)  # catch trailing entry

    results.append((idx, res_time, exp_time))
summary = pd.DataFrame(results, columns=["index", "residence_time_s", "exposure_time_s"])

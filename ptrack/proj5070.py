import sys
import numpy as np
from netCDF4 import Dataset
from pyproj import Transformer

# Usage: python proj5070.py infile.nc outfile.nc

if len(sys.argv) != 3:
    raise SystemExit("Usage: python proj5070.py infile.nc outfile.nc")

infile = sys.argv[1]
outfile = sys.argv[2]

transformer = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)

with Dataset(infile, "r") as src, Dataset(outfile, "w") as dst:

    # Copy dimensions
    for name, dim in src.dimensions.items():
        dst.createDimension(name, None if dim.isunlimited() else len(dim))

    # Copy variables except lon/lat
    for name, var in src.variables.items():
        if name in ("lon", "lat"):
            continue
        out = dst.createVariable(name, var.dtype, var.dimensions, zlib=True)
        out.setncatts(var.__dict__)
        out[:] = var[:]

    # Read lon/lat fully into RAM
    lon = src.variables["lon"][:]
    lat = src.variables["lat"][:]

    # Vectorized projection: runs in compiled C code (very fast)
    x, y = transformer.transform(lon, lat)

    lon_var = dst.createVariable("lon", "f8", src.variables["lon"].dimensions, zlib=True)
    lat_var = dst.createVariable("lat", "f8", src.variables["lat"].dimensions, zlib=True)

    # Add projection metadata
    lon_var.long_name = "x coordinate in EPSG:5070 (NAD83 / Conus Albers)"
    lon_var.units = "meters"
    lat_var.long_name = "y coordinate in EPSG:5070 (NAD83 / Conus Albers)"
    lat_var.units = "meters"
    
    lon_var[:] = x
    lat_var[:] = y

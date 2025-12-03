import sys
from netCDF4 import Dataset
import numpy as np
from pyproj import Transformer

if len(sys.argv) != 3:
    print("Usage: python project_coords.py <input.nc> <output.nc>")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

# Set up transformer (from WGS84 to EPSG:5070)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)

with Dataset(infile, "r") as src:
    # Read original coordinates
    lon = src.variables["lon"][:]
    lat = src.variables["lat"][:]

    # Apply projection
    x_proj, y_proj = transformer.transform(lon[:], lat[:])

    # Copy to new NetCDF
    with Dataset(outfile, "w") as dst:
        # Copy dimensions
        for name, dim in src.dimensions.items():
            dst.createDimension(name, len(dim) if not dim.isunlimited() else None)

        # Copy all variables except lon/lat
        for name, var in src.variables.items():
            if name in ["lon", "lat"]:
                continue
            out_var = dst.createVariable(name, var.datatype, var.dimensions)
            out_var[:] = var[:]
            # Copy attributes
            for attr in var.ncattrs():
                setattr(out_var, attr, getattr(var, attr))

        # Create projected coordinate variables
        x_var = dst.createVariable("x_proj", np.float64, ("pid",))
        y_var = dst.createVariable("y_proj", np.float64, ("pid",))
        x_var[:] = x_proj
        y_var[:] = y_proj
        x_var.units = "meters"
        y_var.units = "meters"

print(f"Projected coordinates written to {outfile}")

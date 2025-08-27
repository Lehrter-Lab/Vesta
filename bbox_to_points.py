import geopandas as gpd
from shapely.geometry import Point
import numpy as np
import pandas as pd
import rasterio

infile  = "bbox_dissolve.zip"
dem     = r"D:\USA\Projects\VESTA\cudem_peconic\mosaicV4.tif"
spacing = 1500  # spacing in meters
outfile = f"grid_{spacing}m.csv"

def generate_points_within_polygon(shapefile, dem, spacing, output_csv):
    # Load polygons
    gdf      = gpd.read_file(shapefile)
    gdf_proj = gdf.to_crs(epsg=5070)

    # Bounds in projected coordinates
    minx, miny, maxx, maxy = gdf_proj.total_bounds

    # Build grid
    xs = np.arange(minx, maxx, spacing)
    ys = np.arange(miny, maxy, spacing)

    points = []
    for x in xs:
        for y in ys:
            p = Point(x + spacing/2, y + spacing/2)
            if gdf_proj.contains(p).any():
                points.append(p)

    # Convert list of Points into GeoDataFrame
    gdf_points = gpd.GeoDataFrame(geometry=points, crs="EPSG:5070")

    with rasterio.open(dem) as src:
            # Reproject to DEM CRS for sampling
            gdf_points_dem = gdf_points.to_crs(src.crs)
            
            # Extract coordinates
            coords = [(pt.x, pt.y) for pt in gdf_points_dem.geometry]
            
            # Sample raster values
            z_vals = [val[0] if val[0] is not None else np.nan 
                      for val in src.sample(coords)]
    
    # Attach Z values
    gdf_points["z"] = z_vals

    # Reproject to NAD83 geographic (EPSG:4269)
    gdf_points = gdf_points.to_crs(epsg=4759)

    # Extract lon/lat and save to CSV
    df_points = pd.DataFrame({"lon": gdf_points.geometry.x,
                              "lat": gdf_points.geometry.y,
                              "z": gdf_points["z"]})
    df_points.to_csv(output_csv, index=False)

    print(f"Generated {len(df_points)} points and saved to {output_csv} (EPSG:4759 lon/lat)")

# Run
generate_points_within_polygon(infile, dem, spacing, outfile)

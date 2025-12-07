import json
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextily as ctx
from pathlib import Path


# --------------------------------------------------------------------
# Convert Julia CSV + meta.json → GeoDataFrame with polygons
# --------------------------------------------------------------------
def df_to_grid_gdf(csv_path, meta_path, crs="EPSG:5070"):
    df = pd.read_csv(csv_path)
    with open(meta_path, "r") as f:
        meta = json.load(f)

    grid = meta["grid_size"]
    min_x, max_x, min_y, max_y = meta["bounds"]

    geoms = []
    for _, r in df.iterrows():
        xi = r.x_bin
        yi = r.y_bin

        # Compute polygon corners
        x0 = min_x + (xi - 1) * grid
        x1 = x0 + grid
        y0 = min_y + (yi - 1) * grid
        y1 = y0 + grid

        poly = Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
        geoms.append(poly)

    gdf = gpd.GeoDataFrame(df, geometry=geoms, crs=crs)
    return gdf


# --------------------------------------------------------------------
# Your heatmap function (unchanged, slightly cleaned)
# --------------------------------------------------------------------
def plot_heatmap(gdf, value_col, title, output_path,
                 crs="EPSG:3857",
                 cmap="viridis_r",
                 log_scale=True,
                 colorbar_label=None,
                 units=3600):

    data = gdf.copy()

    # Convert exposure-like values from seconds → hours
    if ("time" in value_col.lower() or
        "age" in value_col.lower() or
        "exp" in value_col.lower()):
        data[value_col] = data[value_col] / units

    # Remove non-positive if using log scale
    if log_scale:
        data = data[data[value_col] > 0]
        vmin = data[value_col].min()
        vmax = data[value_col].max()
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        vmin = data[value_col].min()
        vmax = data[value_col].max()
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=(12, 8))

    data.to_crs(crs).plot(column=value_col,
                          ax=ax,
                          cmap=cmap,
                          edgecolor="none",
                          norm=norm,
                          legend=False)

    # Add basemap
    ctx.add_basemap(ax, crs=crs, source=ctx.providers.OpenStreetMap.Mapnik)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6)
    cbar.set_label(colorbar_label or title)

    if log_scale:
        ticks = [vmin, np.sqrt(vmin * vmax), vmax]
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels([f"{t:.2f}" for t in ticks])
    else:
        cbar.ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x:.2f}")
        )

    ax.set_axis_off()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


# --------------------------------------------------------------------
# Call
# --------------------------------------------------------------------
# Season
name = "fall"

# Dynamically build names
csv_path   = f"{name}.csv"
meta_path  = f"{name}.meta.json"
mean_name  = f"{name}_mean_exp_time.png"
total_name = f"{name}_total_exp_time.png"

# Build gdf
gdf = df_to_grid_gdf(csv_path, meta_path)

# Plots
plot_heatmap(gdf,
             value_col="mean_exp_time",
             title="Mean Exposure Time (hours)",
             output_path=mean_name,
             log_scale=False,
             colorbar_label="Hours")
plot_heatmap(gdf,
             value_col="total_exp_time",
             title="Total Exposure Time (hours)",
             output_path=total_name,
             log_scale=True,
             colorbar_label="Hours")

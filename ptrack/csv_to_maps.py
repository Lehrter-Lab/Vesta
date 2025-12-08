import json
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextily as ctx
import os


# --------------------------------------------------------------------
# Convert Julia CSV + meta.json → GeoDataFrame with polygons
# --------------------------------------------------------------------
def df_to_grid_gdf(csv_path, meta_path, value_col=None, crs="EPSG:5070"):
    df = pd.read_csv(csv_path)
    with open(meta_path, "r") as f:
        meta = json.load(f)

    grid = meta["grid_size"]
    min_x, max_x, min_y, max_y = meta["bounds"]

    if value_col is not None:
        df = df[~df[value_col].isna()]

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
# Compute seasonal variability using a reference grid
# --------------------------------------------------------------------
def compute_seasonal_variability(csvs, value_col="mean_exp_time", ref_season="fall"):
    """
    Compute standard deviation of a value across multiple seasonal CSVs
    using a reference grid to preserve layout.
    
    Returns a GeoDataFrame with polygons and variability (std_dev) per grid cell.
    """
    # Load reference grid
    ref_gdf = df_to_grid_gdf(csvs[ref_season], f"{ref_season}.meta.json", value_col=None)
    ref_gdf = ref_gdf.set_index(["x_bin", "y_bin"])
    
    seasonal_values = {}
    
    # Load each season and align to reference grid
    for season, csv_path in csvs.items():
        meta_path = f"{season}.meta.json"
        gdf = df_to_grid_gdf(csv_path, meta_path, value_col=value_col)
        gdf_idxed = gdf.set_index(["x_bin", "y_bin"])
        
        # Join seasonal values to reference grid
        seasonal_values[season] = ref_gdf[['geometry']].join(gdf_idxed[[value_col]])[value_col]
    
    # Combine into a DataFrame
    combined = pd.DataFrame(seasonal_values)
    
    # Compute std deviation ignoring NaNs
    ref_gdf['range'] = combined.max(axis=1) - combined.min(axis=1)
    
    # Reset index for plotting
    ref_gdf = ref_gdf.reset_index()
    
    return ref_gdf


# --------------------------------------------------------------------
# Your heatmap function (unchanged, slightly cleaned)
# --------------------------------------------------------------------
def plot_heatmap(gdf, value_col, title, output_path,
                 crs="EPSG:4326",
                 cmap="viridis",
                 log_scale=True,
                 colorbar_label=None,
                 seconds=86400,
                 extent=None):

    data = gdf.copy()

    # Convert exposure-like values from seconds → hours
    if seconds:
        data[value_col] = data[value_col] / seconds
        
    data = data[~data[value_col].isna()]
    data = data[data[value_col] > 0]
    vmin = data[value_col].min()
    vmax = data[value_col].max()
    
    if log_scale:
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=(12, 8))

    data.to_crs(crs).plot(column=value_col,
                          ax=ax,
                          cmap=cmap,
                          edgecolor="black",
                          linewidth=0.1,
                          norm=norm,
                          legend=False)
    
    # Set map frame extent
    if extent:
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
    # Add basemap
    ctx.add_basemap(ax, crs=crs, source=ctx.providers.CartoDB.PositronNoLabels)

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
            plt.FuncFormatter(lambda x, _: f"{x:.2f}"))
    
    # Add ticks
    ax.set_axis_on()
    ax.tick_params(axis='both', which='major', labelsize=10, color='black', length=5)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    
    # Add a black border around plot area
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1)
        
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# --------------------------------------------------------------------
# Call
# --------------------------------------------------------------------
# Ensure the figures folder exists
os.makedirs("figures", exist_ok=True)

# Inputs
seasons = ["fall", "winter", "spring", "summer"]
value_cols = [("mean_time_to_exit", "Mean Residence Time (hours)", "res"),
              ("mean_exp_time", "Mean Exposure Time (hours)", "exp")]

# Get extent for long island
ref = df_to_grid_gdf("fall.csv", "fall.meta.json", value_col=None)
ref = ref.to_crs("EPSG:4326")
xmin, ymin, xmax, ymax = ref.total_bounds
domain = (xmin*1.0005, xmax, ymin, ymax*0.995)
# Get extent for peconic
ref = df_to_grid_gdf("fall.csv", "fall.meta.json", value_col="mean_time_to_exit")
ref = ref.to_crs("EPSG:4326")
xmin, ymin, xmax, ymax = ref.total_bounds
zoom = (xmin*1.0005, xmax*0.9995, ymin*0.9995, ymax*1.0005)

# Loop over seasons and value types
for season in seasons:
    csv_path = f"{season}.csv"
    meta_path = f"{season}.meta.json"
    
    for col, title, suffix in value_cols:
        if suffix == "res":
            extent = zoom
        else:
            extent = domain
        gdf = df_to_grid_gdf(csv_path, meta_path, value_col=col)
        output_path = f"figures/{season}_mean_{suffix}_time.png"
        plot_heatmap(gdf,
                     value_col=col,
                     title=title,
                     output_path=output_path,
                     log_scale=False,
                     colorbar_label="Days",
                     extent=extent)

# Seasonal variability
csvs = {season: f"{season}.csv" for season in seasons}
gdf_var = compute_seasonal_variability(csvs, value_col="mean_exp_time")
plot_heatmap(gdf_var,
             value_col="range",
             title="Seasonal Variability in Mean Exposure Time",
             output_path="figures/seasonal_variability.png",
             log_scale=False,
             colorbar_label="Days",
             extent=extent)

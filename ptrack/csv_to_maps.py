import json
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import contextily as ctx
import os


# --------------------------------------------------------------------
# Convert Julia CSV + meta.json to GeoDataFrame with polygons
# --------------------------------------------------------------------
def df_to_gdf(csv_path, meta_path, value_col=None, crs="EPSG:5070"):
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
    # Determine global bounds and grid size
    bounds_list = []
    grid_sizes = []
    for season, csv_path in csvs.items():
        meta_path = f"{season}.meta.json"
        with open(meta_path, "r") as f:
            meta = json.load(f)
        min_x, max_x, min_y, max_y = meta["bounds"]
        bounds_list.append((min_x, min_y, max_x, max_y))
        grid_sizes.append(meta["grid_size"])
    global_bounds = (min(b[0] for b in bounds_list),
                     min(b[1] for b in bounds_list),
                     max(b[2] for b in bounds_list),
                     max(b[3] for b in bounds_list))
    grid_size = min(grid_sizes)
    
    # Create template grid
    xmin, ymin, xmax, ymax = global_bounds
    x_bins = np.arange(xmin, xmax, grid_size)
    y_bins = np.arange(ymin, ymax, grid_size)
    polygons, x_idx, y_idx = [], [], []
    for i, x0 in enumerate(x_bins):
        for j, y0 in enumerate(y_bins):
            polygons.append(Polygon([(x0, y0), (x0+grid_size, y0),
                                     (x0+grid_size, y0+grid_size), (x0, y0+grid_size)]))
            x_idx.append(i)
            y_idx.append(j)
    template = gpd.GeoDataFrame({'x_idx': x_idx, 'y_idx': y_idx}, 
                                geometry=polygons, crs="EPSG:5070")
    template = template.reset_index(drop=False).rename(columns={'index':'template_idx'})

    
    # Align seasons
    seasonal_arrays = []
    for season, csv_path in csvs.items():
        meta_path = f"{season}.meta.json"
        gdf       = df_to_gdf(csv_path, meta_path, value_col=value_col)
        # Spatial join
        joined = gpd.sjoin(template, gdf[['geometry', value_col]], how='left', predicate='intersects')
        # Aggregate multiple polygons per template cell
        grouped = joined.groupby('template_idx')[value_col].mean()
        seasonal_array = np.full(len(template), np.nan)
        seasonal_array[grouped.index] = grouped.values
        seasonal_arrays.append(seasonal_array)
    
    seasonal_arrays = np.stack(seasonal_arrays, axis=1)

    # Compute range ignoring NaNs
    template['range'] = np.nanmax(seasonal_arrays, axis=1) - np.nanmin(seasonal_arrays, axis=1)
    
    return template


# --------------------------------------------------------------------
# Heatmap function
# --------------------------------------------------------------------
def plot_heatmap(gdf, value_col, output_path,
                 title="",
                 crs="EPSG:4326",
                 cmap="viridis",
                 log_scale=True,
                 colorbar_label=None,
                 seconds=86400,
                 extent=None,
                 vmin=None,
                 vmax=None):

    data = gdf.copy()

    # Convert exposure-like values from seconds → hours
    if seconds:
        data[value_col] = data[value_col] / seconds
        
    data = data[~data[value_col].isna()]
    data = data[data[value_col] > 0]
    
    # Color bar range
    vmin = vmin if vmin is not None else data[value_col].min()
    vmax = vmax if vmax is not None else data[value_col].max()
    if log_scale:
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Plot data on fig
    fig, ax = plt.subplots(figsize=(16, 11))
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
    ctx.add_basemap(ax, crs=crs, 
                    source=ctx.providers.CartoDB.PositronNoLabels, zoom=15)
    
    # Get axes for cbar mapping
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="3%", pad=-1.2)

    # Colorbar
    sm   = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cax)
    
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
        
    # Add season annotation in top-left corner
    ax.text(0.025, 0.97, title, 
            transform=ax.transAxes,   # Use axes coordinates (0 to 1)
            fontsize=20, 
            fontweight='bold',
            va='top',                 # vertical alignment
            ha='left',                # horizontal alignment
            color='black')
    
    fig.tight_layout()
        
    plt.title("")
    plt.tight_layout()
    plt.savefig(output_path, dpi=450)
    plt.close()

# --------------------------------------------------------------------
# Quad Heatmap function
# --------------------------------------------------------------------

def plot_quad_heatmap(seasons, csv_base, value_col, extent, cmap="viridis", 
                      log_scale=False, seconds=86400, vmin=None, vmax=None, output_path="quad.png"):
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for i, season in enumerate(seasons):
        ax = axes[i]
        csv_path = f"{season}.csv"
        meta_path = f"{season}.meta.json"
        gdf = df_to_gdf(csv_path, meta_path, value_col=value_col)
        
        # Scale to hours/days
        if seconds:
            gdf[value_col] = gdf[value_col] / seconds
        
        data = gdf[~gdf[value_col].isna()]
        data = data[data[value_col] > 0]
        
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax) if log_scale else mcolors.Normalize(vmin=vmin, vmax=vmax)
        data.to_crs("EPSG:4326").plot(column=value_col, ax=ax, cmap=cmap,
                                      edgecolor="black", linewidth=0.1, norm=norm, legend=False)
        
        # Set extent
        xmin, xmax, ymin, ymax = extent
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
        # Basemap
        ctx.add_basemap(ax, crs="EPSG:4326", source=ctx.providers.CartoDB.PositronNoLabels)
        
        # Annotation
        ax.text(0.025, 0.97, season.capitalize(),
                transform=ax.transAxes, fontsize=18, fontweight='bold',
                va='top', ha='left', color='black')
        
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(1)
    
    # Add shared colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.02)
    cbar.set_label(f"{value_col} (hours)" if seconds else value_col)
    
    fig.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# --------------------------------------------------------------------
# Call
# --------------------------------------------------------------------
# Ensure the figures folder exists
os.makedirs("figures", exist_ok=True)

# Inputs
seasons = ["winter", "spring", "summer", "fall"]
value_cols = [("mean_time_to_exit", "Mean Residence Time (hours)", "res"),
              ("mean_exp_time", "Mean Exposure Time (hours)", "exp")]

# Get extent for long island
ref = df_to_gdf("fall.csv", "fall.meta.json", value_col=None).to_crs("EPSG:4326")
xmin, ymin, xmax, ymax = ref.total_bounds
domain = (xmin*1.0005, xmax, ymin, ymax*0.995)
# Get extent for peconic
ref = df_to_gdf("fall.csv", "fall.meta.json", value_col="mean_time_to_exit").to_crs("EPSG:4326")
xmin, ymin, xmax, ymax = ref.total_bounds
zoom = (xmin*1.0005, xmax*0.9995, ymin*0.9995, ymax*1.0005)

# Get colorbar range
cbars= {}
for col, _, _ in value_cols:
    all_vals = []
    for season in seasons:
        gdf = df_to_gdf(f"{season}.csv", f"{season}.meta.json", value_col=col)
        all_vals.append(gdf[col].dropna().values)
    all_vals   = np.concatenate(all_vals)
    cbars[col] = (all_vals.min(), all_vals.max())
    print(f"{col} global range: {cbars[col][0]:.2f} : {cbars[col][1]:.2f}")


# Loop over seasons and value types
for season in seasons:
    csv_path  = f"{season}.csv"
    meta_path = f"{season}.meta.json"
    anno      = season.capitalize()
    
    for col, title, suffix in value_cols:
        if suffix == "res":
            extent = zoom
            seconds = 86400
            #cbar    = (None,None)
            cbar    = (cbars[col][0] / seconds, 
                       cbars[col][1] / seconds)
        else:
            extent  = domain
            seconds = 3600
            cbar    = (cbars[col][0] / seconds, 
                       cbars[col][1] / seconds)
        gdf = df_to_gdf(csv_path, meta_path, value_col=col)
        output_path = f"figures/{season}_mean_{suffix}_time.png"
        plot_heatmap(gdf,
                     value_col=col,
                     title=anno,
                     output_path=output_path,
                     log_scale=False,
                     colorbar_label="",
                     seconds=seconds,
                     extent=extent,
                     vmin=cbar[0],
                     vmax=cbar[1])

# Seasonal variability
csvs = {season: f"{season}.csv" for season in seasons}
for col in ["mean_time_to_exit", "mean_exp_time"]:
    gdf_var = compute_seasonal_variability(csvs, value_col=col)
    
    # Set seconds based on type
    if col == "mean_time_to_exit":
        seconds = 86400 
        extent  = zoom 
    else:
        seconds = 3600
        extent  = domain
    
    output_path = f"figures/seasonal_variability_{col}.png"
    plot_heatmap(gdf_var,
                 value_col="range",
                 output_path=output_path,
                 log_scale=False,
                 colorbar_label="",
                 seconds=seconds,
                 extent=extent)
    
# Quad
value_col   = "mean_time_to_exit"
output_path = "figures/quad_mean_res_time.png"

# Use same extent for all subplots
ref = df_to_gdf("fall.csv", "fall.meta.json", value_col=value_col).to_crs("EPSG:4326")
xmin, ymin, xmax, ymax = ref.total_bounds
extent                 = (xmin*1.0005, xmax*0.9995, ymin*0.9995, ymax*1.0005)

# Determine vmin/vmax globally
all_vals = []
for season in seasons:
    gdf = df_to_gdf(f"{season}.csv", f"{season}.meta.json", value_col=value_col)
    all_vals.append(gdf[value_col].dropna().values)
all_vals = np.concatenate(all_vals)
vmin, vmax = all_vals.min()/86400, all_vals.max()/86400  # convert seconds → hours

plot_quad_heatmap(seasons, "base", value_col, extent, cmap="viridis",
                  log_scale=False, seconds=86400, vmin=vmin, vmax=vmax,
                  output_path=output_path)

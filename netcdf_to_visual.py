import xarray as xr
import os
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point, box
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextily as ctx
from pyproj import Geod

# --- Config ---
path_nc     = "particle.nc"
path_domain = "zip:bbox_dissolve.zip"
output_nc   = "particle_enhanced.nc"
output_exp  = "local_exposure_map.png"
output_age  = "local_age_map.png"

def check_enhanced_files(path_enhanced, path_times):
    """Check if enhanced files and variables exist."""
    try:
        if not (os.path.exists(path_enhanced) and os.path.exists(path_times)):
            return False, None, None

        ds_enh = xr.open_dataset(path_enhanced)
        ds_t   = xr.open_dataset(path_times)

        has_inside   = "inside"   in ds_enh.variables
        has_res_time = "res_time" in ds_t.variables
        has_exp_time = "exp_time" in ds_t.variables

        if has_inside and has_res_time and has_exp_time:
            return True, ds_enh, ds_t
        else:
            return False, None, None

    except Exception as e:
        print(f"Error checking enhanced files: {e}")
        return False, None, None

def compute_times(group):
    times      = group["time"].values
    inside     = group["inside"].values.astype(bool)
    first_time = times[0]

    exited_indices = np.where(~inside)[0]
    if exited_indices.size > 0:
        t_exit = times[exited_indices[0]]
    elif inside.any():
        t_exit = times[inside].max()
    else:
        t_exit = np.nan
    res_time = t_exit - first_time if not np.isnan(t_exit) else np.nan

    unique_inside_times = np.unique(times[inside])
    exp_time = np.sum(np.diff(unique_inside_times)) if unique_inside_times.size >= 2 else 0.0

    return np.array([res_time, exp_time])

def compute_local(df, grid_size_deg=0.005, target_crs="EPSG:4326"):
    df = df.copy()
    df['x_bin'] = (df['longitude'] // grid_size_deg).astype(int)
    df['y_bin'] = (df['latitude']  // grid_size_deg).astype(int)

    # Calculate delta t
    df       = df.sort_values(by=['pid', 'time'])
    df['dt'] = df.groupby('pid')['time'].diff().fillna(0)
    
    # Total time each particle spends in each cell including repeat visits
    exposure_per_cell = df.groupby(['x_bin', 'y_bin', 'pid'])['dt'].sum().reset_index()
    
    # Water age weighted by dt per row
    df_valid                  = df[df['dt'] > 0].copy()
    df_valid['time_weighted'] = df_valid['time'] * df_valid['dt']
    
    # Aggregate per cell:
    grouped_exposure = exposure_per_cell.groupby(['x_bin', 'y_bin']).agg({
        'dt': 'mean',
        'pid': 'nunique'}).rename(columns=
                                  {'dt': 'mean_exp_time', 'pid': 'n_particles'})
    grouped_age = df_valid.groupby(['x_bin', 'y_bin']).agg({
        'time_weighted': 'sum',
        'dt': 'sum'})
    grouped = grouped_exposure.join(grouped_age, how='outer').reset_index()
    grouped['mean_water_age'] = grouped['time_weighted'] / grouped['dt']
    
    def bin_to_polygon(xi, yi):
        minx = xi   * grid_size_deg
        maxx = minx + grid_size_deg
        miny = yi   * grid_size_deg
        maxy = miny + grid_size_deg
        return box(minx, miny, maxx, maxy)

    grouped['geometry'] = grouped.apply(lambda row: bin_to_polygon(row['x_bin'], row['y_bin']), axis=1)
    gdf = gpd.GeoDataFrame(grouped, geometry='geometry', crs="EPSG:4326")
    gdf = gdf.to_crs(target_crs)
    return gdf

def compute_mean_distance_per_day(df):
    geod       = Geod(ellps="WGS84")
    df_sorted  = df.sort_values(by=["pid", "time"]).copy()
    df_sorted["lon_prev"] = df_sorted.groupby("pid")["longitude"].shift()
    df_sorted["lat_prev"] = df_sorted.groupby("pid")["latitude"].shift()
    mask_valid = df_sorted["lon_prev"].notna()
    df_valid   = df_sorted[mask_valid].copy()

    _, _, distances = geod.inv(
        df_valid["lon_prev"].values,
        df_valid["lat_prev"].values,
        df_valid["longitude"].values,
        df_valid["latitude"].values)
    df_valid["dxy"]  = distances
    total_distance_m = df_valid.groupby("pid")["dxy"].sum()
    duration_s       = df_sorted.groupby("pid")["time"].agg(lambda x: x.max() - x.min())
    duration_days    = duration_s / 86400
    mean_km_per_day  = (total_distance_m / 1000) / duration_days
    return pd.DataFrame({"mean_km_per_day": mean_km_per_day})

def plot_heatmap(gdf, value_col, title, output_path, crs="EPSG:3857", 
                 cmap="viridis_r", log_scale=True, colorbar_label=None, units=3600):
    # Convert time-like data from seconds to days if relevant
    data = gdf.copy()
    if "time" in value_col or "age" in value_col or "exp" in value_col:
        data[value_col] = data[value_col] / units  # seconds to hours
    
    # Filter out non-positive values if log scale
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
    data.to_crs(crs).plot(
        column=value_col, 
        ax=ax, 
        cmap=cmap, 
        edgecolor='none', 
        norm=norm,
        legend=False)

    ctx.add_basemap(ax, crs=crs, source=ctx.providers.OpenStreetMap.Mapnik)

    # Colorbar
    sm   = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6)
    cbar.set_label(colorbar_label or title)

    if log_scale:
        ticks = [vmin, np.sqrt(vmin * vmax), vmax]
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels([f"{t:.2f}" for t in ticks])
    else:
        cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.2f}"))

    ax.set_axis_off()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

def main():
    # Load domain polygon and unify
    gdf_dom     = gpd.read_file(path_domain).to_crs("EPSG:4326")
    domain_geom = gdf_dom.union_all()
    # Load dataset with chunking for memory efficiency
    ds = xr.open_dataset(path_nc, chunks={"pid": 1000})
    # Check if enhanced files and variables exist
    found_all, ds_enhanced, ds_times = check_enhanced_files(output_nc, "times.nc")

    if found_all:
        print("Enhanced variables found — skipping computation.")
        ds        = ds_enhanced
        df        = ds.to_dataframe().reset_index()
        times_arr = pd.DataFrame({
            "res_time": ds_times["res_time"].values,
            "exp_time": ds_times["exp_time"].values
        }, index=ds_times.coords["pid"].values)

    else:
        print("Enhanced variables not found — computing inside, residence, and exposure times.")
        df = ds.to_dataframe().reset_index()
        # Filter out static particles
        moved = df.groupby("pid")[["longitude", "latitude"]].transform(lambda x: x.nunique() > 1)
        df    = df[(moved["longitude"]) | (moved["latitude"])].copy()
        # Bounding box of the domain (minx, miny, maxx, maxy)
        minx, miny, maxx, maxy = domain_geom.bounds
        # Pre-filter: only keep points inside bounding box
        bbox_mask = ((df["longitude"] >= minx) & (df["longitude"] <= maxx) &
                     (df["latitude"]  >= miny) & (df["latitude"]  <= maxy))
        # Initialize inside flags to False, and check only for candidates in bbox
        inside_flags = np.zeros(len(df), dtype=bool)
        # Check contains() only for candidates
        candidate_idx = np.where(bbox_mask)[0]
        points_subset = [Point(xy) for xy in zip(df.loc[candidate_idx, "longitude"],
                                                 df.loc[candidate_idx, "latitude"])]
        inside_checked = np.array([domain_geom.contains(pt) for pt in points_subset])
        # Assign results
        inside_flags[candidate_idx] = inside_checked
        df["inside"]                = inside_flags
        # Compute residence and exposure times
        times_arr = df.groupby("pid", observed=True, group_keys=False).apply(compute_times)
        times_arr = pd.DataFrame(times_arr.tolist(), index=times_arr.index, columns=["res_time", "exp_time"])
        # Assign back the new variables to xarray Dataset and save enhanced file
        da_inside = xr.DataArray(inside_flags, dims=["point"])
        ds        = ds.assign({"inside": da_inside})
        ds.to_netcdf(output_nc)
        # Assign to a new times.nc
        da_res = xr.DataArray(
            times_arr['res_time'].values,
            coords={"pid": times_arr.index.values},  # 9600 particle IDs
            dims=["pid"],
            name="res_time")
        da_exp = xr.DataArray(
            times_arr['exp_time'].values,
            coords={"pid": times_arr.index.values},
            dims=["pid"],
            name="exp_time")
        ds_times = xr.Dataset({
            "res_time": da_res,
            "exp_time": da_exp})
        ds_times.to_netcdf("times.nc")
    # Print average stats for filtered particles
    print(f"Average residence time: {times_arr['res_time'].mean()/86400:.2f} days")
    print(f"Average exposure time: {times_arr['exp_time'].mean()/86400:.2f} days")
    # Merge exposure times back for binning and plotting
    df  = df.merge(times_arr[['res_time', 'exp_time']], left_on='pid', right_index=True, how='left')
    # Compute daily particle transit
    df_speed = compute_mean_distance_per_day(df)
    print(f"Average travel speed: {df_speed['mean_km_per_day'].mean():.2f} km/day")
    # Compute local exposure and plot
    gdf = compute_local(df, grid_size_deg=0.005)
    plot_heatmap(gdf, value_col='mean_exp_time', output_path=output_exp,
                   title="Local Exposure Time (hours)", log_scale=True)
    # Compute local water age and plot
    plot_heatmap(gdf, value_col='mean_water_age', output_path=output_age,
                   title="Local Water Age (days)", log_scale=False, units=86400)
    # Save centroids
    # Project to EPSG:5070 to compute accurate centroids
    gdf_proj       = gdf.to_crs("EPSG:5070")
    centroids_proj = gdf_proj.geometry.centroid
    # Convert centroids back to EPSG:4326
    centroids_geo = gpd.GeoSeries(centroids_proj, crs="EPSG:5070").to_crs("EPSG:4326")
    # Create a GeoDataFrame of centroids with original attributes
    gdf_centroids = gpd.GeoDataFrame(gdf.copy(), geometry=centroids_geo, crs="EPSG:4326")
    # Mask: keep only centroids inside the domain polygon
    mask_inside = gdf_centroids.geometry.within(domain_geom)
    gdf_masked  = gdf_centroids[mask_inside].copy()
    # Extract x/y and save to CSV
    gdf_masked["x"] = gdf_masked.geometry.x
    gdf_masked["y"] = gdf_masked.geometry.y
    gdf_masked[["x", "y"]].to_csv("grid_centers_xy.csv", index=False)
    
if __name__ == "__main__":
    main()

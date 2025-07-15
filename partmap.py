# CMikolaitis @ USA/DISL

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
import contextily as ctx
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyproj

# Input files
infile     = 'particle.pth'
bboxfile   = 'zip:bbox_dissolve.zip'
outmap     = 'particle_times.pdf'
outtrack   = 'particle_track.pdf'

# Plotting parameters
grid_res       = 0.005
figsize        = (16, 8)
cmap           = 'plasma'
cbar_label     = 'Time (days)'
titles         = ['Mean Residence Time (days)', 'Mean Exposure Time (days)']
domain_lw      = 1
domain_ec      = 'black'
cbar_width     = "5%"
cbar_pad       = 0.05
tick_decimals  = 2
font_size      = 12
basemap_source = ctx.providers.OpenStreetMap.Mapnik
dpi            = 300

# Transformer from EPSG:3857 to EPSG:4326 for tick labels
projector = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

# Read particle.pth format
def read_particles(filename):
    with open(filename) as f:
        lines = f.readlines()
    times, lons, lats, ids = [], [], [], []
    i = 0
    while i < len(lines):
        parts = lines[i].strip().split()
        if len(parts) == 2:
            t, n = float(parts[0]), int(parts[1])
            for j in range(n):
                pid, lon, lat = lines[i + 1 + j].strip().split()[:3]
                times.append(t)
                ids.append(int(pid))
                lons.append(float(lon))
                lats.append(float(lat))
            i += n + 1
        else:
            i += 1
    return np.array(lons), np.array(lats), np.array(times), np.array(ids)

# Prepare residence and exposure times combined in one GeoDataFrame
def prepare_data(lons, lats, times, ids, bbox, skip_static=False):
    df_raw = pd.DataFrame({'id': ids, 'lon': lons, 'lat': lats, 'time': times})
    
    # If particle never moves, omit
    if skip_static:
        # Find particles with no movement (all lon and lat identical)
        grouped = df_raw.groupby('id')
        # Compute range of lon and lat per id
        lon_range = grouped['lon'].transform(lambda x: x.max() - x.min())
        lat_range = grouped['lat'].transform(lambda x: x.max() - x.min())
        # Keep only moving particles
        moving_mask = (lon_range > 0) | (lat_range > 0)
        df_raw = df_raw[moving_mask].copy()
        
    # Convert to geometry 
    geom = gpd.points_from_xy(df_raw.lon, df_raw.lat)
    df_geo = gpd.GeoDataFrame(df_raw, geometry=geom, crs=bbox.crs)
    domain_geom = bbox.union_all()
    df_geo['inside'] = df_geo.geometry.within(domain_geom)
    
    # Residence time calculation:
    # For each id, find min time outside domain (if any), else max time inside domain.
    # Compute first appearance time for each id to subtract later.
    # Get min time outside domain per id (NaN if no outside time)
    outside_times = (
        df_geo.loc[~df_geo['inside']]
              .groupby('id')['time']
              .min()
              .rename('min_outside_time')
    )
    # Get max time inside domain per id
    inside_times = (
        df_geo.loc[df_geo['inside']]
              .groupby('id')['time']
              .max()
              .rename('max_inside_time')
    )
    # Get first appearance time per id
    first_times = df_geo.groupby('id')['time'].min().rename('first_time')
    # Combine into one DataFrame
    df_res = pd.DataFrame({
        'residence_time': outside_times.combine_first(inside_times),
    })
    df_res['residence_time'] = df_res['residence_time'] - first_times
    df_res = df_res.reset_index()
    
    # Exposure time calculation:
    # Total time particle spends inside domain, computed as sum of time diffs in sorted unique times inside domain per id
    def exposure_time(times_array):
        if len(times_array) <= 1:
            return 0
        sorted_t = np.sort(times_array)
        return np.sum(np.diff(sorted_t))
    exposure_times = (
        df_geo[df_geo['inside']]
        .groupby('id')['time']
        .agg(lambda x: exposure_time(x.unique()))
        .rename('exposure_time')
        .reset_index()
    )
    # Injection site (first lon/lat per id)
    df_rel = df_raw.groupby('id')[['lon', 'lat']].first().reset_index()
    
    # Merge all data
    df_out = df_res.merge(exposure_times, on='id', how='left')
    df_out = df_out.merge(df_rel, on='id', how='left')
    df_out['exposure_time'] = df_out['exposure_time'].fillna(0)
    # GeoDataFrame with injection site geometry
    df_out = gpd.GeoDataFrame(df_out, geometry=gpd.points_from_xy(df_out.lon, df_out.lat), crs=bbox.crs)
    return df_raw, df_out

# Bin and project to EPSG:3857
def bin_and_project(df, bbox, value_col, res=None, method='mean', bins=None):
    if res is not None:
        minx, miny, maxx, maxy = bbox.total_bounds
        lon_bins = np.arange(minx, maxx + res, res)
        lat_bins = np.arange(miny, maxy + res, res)
        bins = [lon_bins, lat_bins]
    stat, x_edges, y_edges, _ = binned_statistic_2d(
        df['lon'], df['lat'], df[value_col], statistic=method, bins=bins)
    xs = (x_edges[:-1] + x_edges[1:]) / 2
    ys = (y_edges[:-1] + y_edges[1:]) / 2
    xv, yv = np.meshgrid(xs, ys)
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xv.ravel(), yv.ravel()), crs=bbox.crs)
    proj = pts.to_crs(epsg=3857)
    return stat, proj.geometry.x.values.reshape(xv.shape), proj.geometry.y.values.reshape(yv.shape)

# Plot residence + exposure heatmaps side-by-side
def plot_heatmaps(stat_res, stat_exp, x, y, bbox):
    plt.rcParams.update({'font.size': font_size})
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    # Convert seconds to days
    stat_res_days = stat_res / 86400
    stat_exp_days = stat_exp / 86400
    vmin, vmax = np.nanmin([stat_res_days, stat_exp_days]), np.nanmax([stat_res_days, stat_exp_days])

    for ax, stat, title in zip(axes, [stat_res_days, stat_exp_days], titles):
        im = ax.pcolormesh(x, y, stat.T, cmap=cmap, vmin=vmin, vmax=vmax)
        bbox.boundary.plot(ax=ax, edgecolor=domain_ec, linewidth=domain_lw)
        ctx.add_basemap(ax, crs='EPSG:3857', source=basemap_source)
        ax.set_title(title)
        ax.tick_params(axis='y', labelleft=True)
        xticks = ax.get_xticks()
        yticks = ax.get_yticks()
        lon_ticks, _ = projector.transform(xticks, np.zeros_like(xticks))
        _, lat_ticks = projector.transform(np.zeros_like(yticks), yticks)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels([f"{lon:.{tick_decimals}f}" for lon in lon_ticks])
        ax.set_yticklabels([f"{lat:.{tick_decimals}f}" for lat in lat_ticks])

        cax = make_axes_locatable(ax).append_axes("right", size=cbar_width, pad=cbar_pad)
        fig.colorbar(im, cax=cax, label=cbar_label)

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    fig.savefig(outmap, dpi=dpi, bbox_inches='tight')
    plt.show()

# Plot heatmap of time spent per cell for all paths (residence) with lat/lon tick labels
def plot_paths_heatmap(lons, lats, times, ids, bbox, log=False):
    df_all = pd.DataFrame({'id': ids, 'lon': lons, 'lat': lats, 'time': times})
    df_all = df_all.sort_values(['id', 'time'])
    df_all['dt'] = df_all.groupby('id')['time'].diff().fillna(0)

    # Bins over all particle data extent
    lon_bins = np.arange(df_all['lon'].min(), df_all['lon'].max() + grid_res, grid_res)
    lat_bins = np.arange(df_all['lat'].min(), df_all['lat'].max() + grid_res, grid_res)

    df_all['lon_bin'] = np.digitize(df_all['lon'], lon_bins) - 1
    df_all['lat_bin'] = np.digitize(df_all['lat'], lat_bins) - 1

    valid = (
        (df_all['lon_bin'] >= 0) & (df_all['lat_bin'] >= 0) &
        (df_all['lon_bin'] < len(lon_bins) - 1) & (df_all['lat_bin'] < len(lat_bins) - 1)
    )
    df_all = df_all[valid]

    df_bins = df_all.groupby(['lon_bin', 'lat_bin'])['dt'].sum().reset_index(name='dt_sum')

    df_bins['lon'] = lon_bins[df_bins['lon_bin']] + grid_res / 2
    df_bins['lat'] = lat_bins[df_bins['lat_bin']] + grid_res / 2

    stat, x_proj, y_proj = bin_and_project(df_bins, bbox, 'dt_sum', bins=[lon_bins, lat_bins])

    stat_days = stat / 86400
    x_min, x_max = x_proj.min(), x_proj.max()
    y_min, y_max = y_proj.min(), y_proj.max()

    lo, hi = np.nanpercentile(stat_days, [0.5, 99.5])

    if log:
        positive = stat_days[stat_days > 0]
        norm = mpl.colors.LogNorm(vmin=positive.min(), vmax=positive.max(), clip=True)
        cbar_label_local = "Time per grid cell (days, log scale)"
        title = "Time spent per grid cell (days, log scale)"
    else:
        norm = mpl.colors.Normalize(vmin=lo, vmax=hi, clip=True)
        cbar_label_local = "Time spent per grid cell (days)"
        title = "Time spent per grid cell (days)"

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.pcolormesh(x_proj, y_proj, stat_days.T, cmap=cmap, norm=norm)

    bbox.to_crs(epsg=3857).boundary.plot(ax=ax, edgecolor=domain_ec, linewidth=domain_lw)
    ctx.add_basemap(ax, crs='EPSG:3857', source=basemap_source, reset_extent=False)
    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    ax.set_title(title)
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    lon_ticks, _ = projector.transform(xticks, np.zeros_like(xticks))
    _, lat_ticks = projector.transform(np.zeros_like(yticks), yticks)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{lon:.{tick_decimals}f}" for lon in lon_ticks])
    ax.set_yticklabels([f"{lat:.{tick_decimals}f}" for lat in lat_ticks])
    
    cax = make_axes_locatable(ax).append_axes("right", size=cbar_width, pad=cbar_pad)
    fig.colorbar(im, cax=cax, label=cbar_label_local)
    
    fig.savefig(outtrack, dpi=dpi, bbox_inches='tight')
    plt.show()
    
## Main run -------------------------------------------------------------------
bbox                     = gpd.read_file(bboxfile)
lons, lats, times, ids   = read_particles(infile)
df_raw, df_times         = prepare_data(lons, lats, times, ids, bbox)
stat_res, x_proj, y_proj = bin_and_project(df_times, bbox, 'residence_time', res=grid_res)
stat_exp, _, _           = bin_and_project(df_times, bbox, 'exposure_time', res=grid_res)
bbox_proj                = bbox.to_crs(epsg=3857)
# Plot
plot_heatmaps(stat_res, stat_exp, x_proj, y_proj, bbox_proj)
plot_paths_heatmap(lons, lats, times, ids, bbox, True)

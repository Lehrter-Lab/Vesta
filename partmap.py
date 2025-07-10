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

# Prepare residence/exposure data and filtered GeoDataFrame
def prepare_data(lons, lats, times, ids, bbox):
    df_raw = pd.DataFrame({'id': ids, 'lon': lons, 'lat': lats, 'time': times})
    geom   = gpd.points_from_xy(df_raw.lon, df_raw.lat)
    df_geo = gpd.GeoDataFrame(df_raw, geometry=geom, crs=bbox.crs)
    df_geo['inside'] = df_geo.geometry.within(bbox.union_all())
    
    # Time until particle first leaves domain (seconds)
    df_res = df_geo.sort_values('time').groupby('id', group_keys=False).apply(
        lambda g: pd.Series({'residence_time': g.loc[~g.inside, 'time'].min()
                             if (~g.inside).any() else g.time.max()}),
        include_groups=False).reset_index()
    
    # Total time particle spends inside domain (seconds)
    df_exp = df_geo[df_geo.inside].groupby('id')['time'].unique().apply(
        lambda t: np.diff(np.sort(t)).sum() if len(t) > 1 else 0
    ).reset_index(name='exposure_time')
    
    # Get injection site and add to residence and exposure dfs
    df_rel = df_raw.groupby('id')[['lon', 'lat']].first().reset_index()
    df_res = df_res.merge(df_rel, on='id')
    df_exp = df_exp.merge(df_rel, on='id')

    return df_raw, df_res, df_exp

# Bin and project to EPSG:3857
def bin_and_project(df, bbox, value_col, res=None, method='mean', bins=None):
    if res is not None:
        minx, miny, maxx, maxy = bbox.total_bounds
        lon_bins               = np.arange(minx, maxx + res, res)
        lat_bins               = np.arange(miny, maxy + res, res)
        bins                   = [lon_bins, lat_bins]
    stat, x_edges, y_edges, _  = binned_statistic_2d(
        df['lon'], df['lat'], df[value_col], statistic=method, bins=bins)
    xs     = (x_edges[:-1] + x_edges[1:]) / 2
    ys     = (y_edges[:-1] + y_edges[1:]) / 2
    xv, yv = np.meshgrid(xs, ys)
    pts    = gpd.GeoDataFrame(geometry=gpd.points_from_xy(xv.ravel(), yv.ravel()), crs=bbox.crs)
    proj   = pts.to_crs(epsg=3857)
    return stat, proj.geometry.x.values.reshape(xv.shape), proj.geometry.y.values.reshape(yv.shape)

# Plot residence + exposure heatmaps side-by-side
def plot_heatmaps(stat_res, stat_exp, x, y, bbox):
    plt.rcParams.update({'font.size': font_size})
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    # Convert seconds to days for display
    stat_res_days = stat_res / 86400
    stat_exp_days = stat_exp / 86400
    vmin, vmax = np.nanmin([stat_res_days, stat_exp_days]), np.nanmax([stat_res_days, stat_exp_days])

    for ax, stat, title in zip(axes, [stat_res_days, stat_exp_days], titles):
        im = ax.pcolormesh(x, y, stat.T, cmap=cmap, vmin=vmin, vmax=vmax)
        bbox.boundary.plot(ax=ax, edgecolor=domain_ec, linewidth=domain_lw)
        ctx.add_basemap(ax, crs='EPSG:3857', source=basemap_source)
        ax.set_title(title)
        ax.tick_params(axis='y', labelleft=True)
        # Set tick labels to lat/lon (decimal degrees)
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
def plot_paths_heatmap(lons, lats, times, ids, bbox):
    df_all = pd.DataFrame({'id': ids, 'lon': lons, 'lat': lats, 'time': times})
    df_all = df_all.sort_values(['id', 'time'])
    df_all['dt'] = df_all.groupby('id')['time'].diff().fillna(0)

    # Bins over all particle data extent
    lon_bins = np.arange(df_all['lon'].min(), df_all['lon'].max() + grid_res, grid_res)
    lat_bins = np.arange(df_all['lat'].min(), df_all['lat'].max() + grid_res, grid_res)
    # Digitize lon/lat into bins
    df_all['lon_bin'] = np.digitize(df_all['lon'], lon_bins) - 1
    df_all['lat_bin'] = np.digitize(df_all['lat'], lat_bins) - 1
    # Filter valid bins
    valid = (df_all['lon_bin'] >= 0) & (df_all['lat_bin'] >= 0) & \
            (df_all['lon_bin'] < len(lon_bins)-1) & (df_all['lat_bin'] < len(lat_bins)-1)
    df_all = df_all[valid]
    # Sum time spent per bin
    df_bins = df_all.groupby(['lon_bin', 'lat_bin'])['dt'].sum().reset_index(name='dt_sum')
    # Convert bin indices back to lon/lat centers
    df_bins['lon'] = lon_bins[df_bins['lon_bin']] + grid_res/2
    df_bins['lat'] = lat_bins[df_bins['lat_bin']] + grid_res/2
    # Use bin_and_project to get projected coordinates and heatmap
    stat, x_proj, y_proj = bin_and_project(df_bins, bbox, 'dt_sum', bins=[lon_bins, lat_bins])
    # Convert seconds to days
    stat_days = stat / 86400
    # Match figure box to project extent
    x_min, x_max = x_proj.min(), x_proj.max()
    y_min, y_max = y_proj.min(), y_proj.max()
    # Get confidence interval
    lo, hi = np.nanpercentile(stat_days, [0.5, 99.5])
    norm   = mpl.colors.Normalize(vmin=lo, vmax=hi, clip=True)
    # Actually plot
    fig, ax = plt.subplots(figsize=(10, 10))
    im      = ax.pcolormesh(x_proj, y_proj, stat_days.T, cmap=cmap, norm=norm)
    bbox.to_crs(epsg=3857).boundary.plot(ax=ax, edgecolor=domain_ec, linewidth=domain_lw)
    ctx.add_basemap(ax, crs='EPSG:3857', source=basemap_source, reset_extent=False)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    # Titles and ticks
    ax.set_title("Time spent per grid cell (days)")
    xticks       = ax.get_xticks()
    yticks       = ax.get_yticks()
    lon_ticks, _ = projector.transform(xticks, np.zeros_like(xticks))
    _, lat_ticks = projector.transform(np.zeros_like(yticks), yticks)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{lon:.{tick_decimals}f}" for lon in lon_ticks])
    ax.set_yticklabels([f"{lat:.{tick_decimals}f}" for lat in lat_ticks])

    cax = make_axes_locatable(ax).append_axes("right", size=cbar_width, pad=cbar_pad)
    fig.colorbar(im, cax=cax, label=cbar_label)

    fig.savefig(outtrack, dpi=dpi, bbox_inches='tight')
    plt.show()

# Run everything
bbox                     = gpd.read_file(bboxfile)
lons, lats, times, ids   = read_particles(infile)
df_raw, df_res, df_exp   = prepare_data(lons, lats, times, ids, bbox)
stat_res, x_proj, y_proj = bin_and_project(df_res, bbox, 'residence_time', res=grid_res)
stat_exp, _, _           = bin_and_project(df_exp, bbox, 'exposure_time', res=grid_res)
bbox_proj                = bbox.to_crs(epsg=3857)

plot_heatmaps(stat_res, stat_exp, x_proj, y_proj, bbox_proj)
plot_paths_heatmap(lons, lats, times, ids, bbox)

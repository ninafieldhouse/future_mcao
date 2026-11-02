# %%
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from lagranto.plotting import plot_trajs
from lagranto import Tra 
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
from matplotlib import cm
import cmcrameri.cm as cmc
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
import pandas as pd
import netCDF4 as nc
import numpy as np
import datetime as dt
import xarray as xr
from shapely.geometry import Point
from shapely.prepared import prep
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
import config

# %%

# plt.rcParams.update({
#     'font.size': 16,           # Default font size for ticks
#     'axes.titlesize': 20,      # Title font size
#     'axes.labelsize': 20,      # X and Y label font size
#     'legend.fontsize': 20,     # Legend font size
#     'xtick.labelsize': 16,     # X-axis tick font size
#     'ytick.labelsize': 16      # Y-axis tick font size
# })

# plt.rcParams.update({
#     'font.size': 8,           # Default font size for ticks
#     'axes.titlesize': 8,      # Title font size
#     'axes.labelsize': 8,      # X and Y label font size
#     'legend.fontsize': 6,     # Legend font size
#     'xtick.labelsize': 6,     # X-axis tick font size
#     'ytick.labelsize': 6      # Y-axis tick font size
# })

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
# Configuration
map_extent = {
    "JS": [60, 260, 25, 90],  # for map
    "NA": [130, 330, 25, 90]  # for map, North Atlantic in 0–360°E
}

# map_extent = {
#     "JS": [60, 210, 30, 90],  # Tightened from [60, 260, 25, 70]
#     "NA": [160, 320, 40, 90]  # Tightened from [130, 330, 25, 90]
# }

map_extent_DI = {
    "JS": [60, 260, 25, 90],  # for map
    "NA": [130, 330, 25, 90]  # for map, North Atlantic in 0–360°E
}
box = {
    "JS": [130, 150, 30, 60],  # for mask
    "NA": [285, 295, 30, 45]   # for mask
}
central_longitudes = {
    "JS": 120,  # for map projection
    "NA": 270   # for map projection
}
# central_longitudes = {
#     "JS": 135,  # Center of 90-180° range
#     "NA": 240   # Center of 180-300° range (or -120° in -180 to 180 system)
# }
central_longitudes_DI = {
    "JS": 120,  # for map projection
    "NA": 270   # for map projection
}
selected_trajs = {
    "JS": 1000,
    "NA": 100
}
units_dict = {
    'p': '[hPa]',
    'TH': '[K]',
    'SST': '[K]',
    'dTH': '[K 6h$^{-1}$]',
    'dTH_24': '[K 24h$^{-1}$]',
    'COUNT': '[# trajectories]',
    'RESIDENCE': '[hours]'
}

cmaps_dict = {
    'p': plt.get_cmap('Spectral'),
    'TH': plt.get_cmap('Spectral_r'),
    'SST': plt.get_cmap('Spectral_r'),
    'dTH': plt.get_cmap('RdBu_r'),
    'dTH_24': plt.get_cmap('RdBu_r'),
    'COUNT': plt.get_cmap('plasma'),
    'RESIDENCE': plt.get_cmap('viridis')
}

levels_dict = {
    'p': np.arange(400, 1051, 50),
    'TH': np.arange(230, 301, 5),
    'SST': np.arange(260, 301, 5),
    'dTH': np.arange(-1, 1.1, 0.1),
    'dTH_24': np.arange(-3, 3.25, 0.25),
    'COUNT': np.arange(0, 601, 100),  # logarithmic levels for counts
    'RESIDENCE': np.arange(0, 25, 1)
}

path = f'{config.DENSITY_DIR}/'
fig_path = f'{config.PLOTS_ROOT}/density/'

regions = ['JS', 'NA']
timeframes = ['present', 'future']
airstreams = ['ALL', 'BCAR', 'DI', 'nonBCAR', 'LL']
plotvars = ['p', 'TH', 'dTH', 'COUNT', 'RESIDENCE']

# %%
def load_density_data(density_file):
    """Load density netCDF file and return xarray dataset"""
    try:
        ds = xr.open_dataset(density_file, decode_times=True)
        print(f"Loaded density file: {density_file}")
        print(f"Available variables: {list(ds.data_vars.keys())}")
        print(f"Available times: {ds.time.values if 'time' in ds.dims else 'No time dimension'}")
        return ds
    except Exception as e:
        print(f"Error loading {density_file}: {e}")
        return None
    
def mask_data(ds, varname, time_idx):
    data = ds.isel(time=time_idx)[varname]
    data_masked = np.where(data != 0, data, np.nan)  # Mask out zero values
    
    return data_masked

def create_map_projection(region, central_longitudes):
    """Create map projection for given region"""
    map_projection = ccrs.AzimuthalEquidistant(central_longitude=central_longitudes[region], central_latitude=75)
    return map_projection

def create_map_projection(region, central_longitudes):
    """Create map projection for given region"""
    map_projection = ccrs.LambertAzimuthalEqualArea(central_longitude=central_longitudes[region], central_latitude=60)
    return map_projection

def setup_map_axes(ax, region, projection, map_extent):
    """Setup map axes with coastlines, borders, etc."""

    ax.add_feature(cfeature.LAND, facecolor=(217/256, 217/256, 217/256))
    # ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    # ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    ax.set_extent(map_extent[region], crs=ccrs.PlateCarree())  # extent as box: [115,180,20,60]
    ax.set_facecolor("white")

    gl = ax.gridlines(draw_labels=True, dms=True)  # add gridlines for latitude and longitude
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True

    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
    gl.ylocator = mticker.FixedLocator(np.arange(0, 100, 20))

    gl.xlabel_style = {"size": 20, "color": "k"}
    gl.ylabel_style = {"size": 20, "color": "k"}
    
    return ax

def get_colorbar_shrink(extent, base_width=12.0, base_shrink=0.6):
    """
    Calculates colorbar shrink based on a linear scaling formula.

    Parameters:
        fig_width: float
            Width of the figure (in inches)
        base_width: float
            Width at which shrink = base_shrink
        base_shrink: float
            Shrink at the base width

    Returns:
        float: colorbar shrink factor, clipped to a maximum of 1.0 and rounded to 2 decimals
    """
    fig_width = get_figsize_from_map_extent(extent)[0]
    shrink = base_shrink * (base_width / fig_width)
    return min(1.0, round(shrink, 2))


def timedelta_to_hours(td):
    """Convert numpy.timedelta64 to hours as float"""
    return td / np.timedelta64(1, 'h')

def count_data(ds):
    # COUNTS
    counts = mask_data(ds, 'COUNT', time_idx)
    # Use count_data directly (it's already a numpy array)
    count_data = counts
    lon = ds.longitude
    lat = ds.latitude
    total_count = count_data

    # Calculate percentage of total counts for each grid point
    total_counts = np.nansum(total_count)  # Sum all counts, ignoring NaN values
    percentage_grid = (total_count / total_counts) * 100

    # Create a cumulative percentage grid by sorting grid points by their percentage
    # and calculating cumulative sum
    flat_percentages = percentage_grid.flatten()
    flat_percentages = flat_percentages[~np.isnan(flat_percentages)]  # Remove NaN values
    sorted_percentages = np.sort(flat_percentages)[::-1]  # Sort descending

    # Calculate cumulative sum of percentages
    cumsum_percentages = np.cumsum(sorted_percentages)

    # Find threshold values for 50%, 90%, and 95% of total trajectories
    threshold_50 = sorted_percentages[np.where(cumsum_percentages >= 50)[0][0]]
    threshold_90 = sorted_percentages[np.where(cumsum_percentages >= 90)[0][0]]
    threshold_95 = sorted_percentages[np.where(cumsum_percentages >= 95)[0][0]]

    return percentage_grid, threshold_50, threshold_90, threshold_95

def var_data(ds, varname):
    """Extract variable data from dataset"""
    data = mask_data(ds, varname)
    lon_2d, lat_2d = np.meshgrid(ds.longitude, ds.latitude)
    lon_flat = lon_2d.flatten()
    lat_flat = lat_2d.flatten()
    
    # Flatten the arrays for scatter plot
    data_flat = data.flatten()
    
    # Remove NaN values
    mask = ~np.isnan(data_flat)
    lon_points = lon_flat[mask]
    lat_points = lat_flat[mask]
    data_points = data_flat[mask]
    
    return lon_points, lat_points, data_points

def get_figsize_from_map_extent(map_extent, base_height=10):
    """
    Returns a recommended figsize (width, height) for a Cartopy map
    based on a 4-element extent list: [lon_min, lon_max, lat_min, lat_max].
    Handles dateline crossing and rounds width to the nearest inch.

    Parameters:
        map_extent: list or tuple of [lon_min, lon_max, lat_min, lat_max]
        base_height: float (default = 4)
            Desired figure height in inches

    Returns:
        (width, height): tuple of floats
            Width and height for matplotlib figsize
    """
    if len(map_extent) != 4:
        raise ValueError("map_extent must be a list or tuple of [lon_min, lon_max, lat_min, lat_max]")

    lon_min, lon_max, lat_min, lat_max = map_extent

    # Handle dateline crossing
    if lon_max < lon_min:
        lon_span = (360 + lon_max) - lon_min
    else:
        lon_span = lon_max - lon_min

    lat_span = lat_max - lat_min
    if lat_span == 0:
        raise ValueError("Latitude span cannot be zero.")

    aspect_ratio = lon_span / lat_span
    raw_width = aspect_ratio * base_height

    width = round(raw_width)  # round to nearest inch
    return (width, base_height)
#%% Plot present and future
regions = ['JS', 'NA']
regions = ['JS']
timeframes = ['present', 'future']
# airstreams = ['BCAR', 'DI', 'nonBCAR', 'LL', 'ALL']
airstreams = ['ALL', 'UL', 'LL']  # For now, only plot ALL airstream
airstreams = ['LL']  # For now, only plot ALL airstreams
# airstreams = ['ALL', 'BCAR', 'DI']  # For now, only plot BCAR airstream

#density files:
density_types = ['all']

# density_types = ['timesteps', 'all']
for density_type in density_types:
    if density_type == 'all':
        levels_dict['dTH_24'] = np.arange(-1, 1.1, 0.1)  # Adjust levels for 24-hourly dTH
    else:
        levels_dict['dTH_24'] = np.arange(-3, 3.25, 0.25)
    for region in regions:
        for timeframe in timeframes:
            for airstream in airstreams:
                density_file = f"{path}density_{density_type}_CAO_{airstream}_trajs_{region}_{timeframe}_ALL_with_dTH"
                print(f"Processing density file: {density_file}")
                ds = load_density_data(density_file)
                lon = ds.longitude
                lat = ds.latitude

                if 'time' in ds.dims:
                    time_hours = ds.time.values / np.timedelta64(1, 'h')
                    print("Time coordinate in hours:", time_hours)
                
                # Create meshgrid for coordinates
                lon_2d, lat_2d = np.meshgrid(ds.longitude, ds.latitude)
                lon_flat = lon_2d.flatten()
                lat_flat = lat_2d.flatten()

                varname = 'dTH_24'
                cmap = cmaps_dict[varname]
                levels = levels_dict[varname]
                time_idxs = np.arange(len(ds.time)) if 'time' in ds.dims else [0]
                for time_idx in time_idxs:
                    # DATA
                    data = mask_data(ds, varname, time_idx)
                    #COUNTS
                    percentage_grid, threshold_50, threshold_90, threshold_95 = count_data(ds)

                    mask_95 = np.where(percentage_grid<=threshold_95, np.nan, 1)

                    # Apply mask to data
                    data = data * mask_95

                    # Flatten the arrays for scatter plot
                    data_flat = data.flatten()

                    # Remove NaN values
                    mask = ~np.isnan(data_flat)
                    lon_points = lon_flat[mask]
                    lat_points = lat_flat[mask]
                    data_points = data_flat[mask]

                    if airstream in ['ALL', 'UL', 'LL']:
                        map_extent = map_extent
                        central_longitudes = central_longitudes
                    else:
                        map_extent = map_extent_DI
                        central_longitudes = central_longitudes_DI
                    fig = plt.figure(figsize=(10, 8), layout='constrained')
                    fig = plt.figure(figsize=get_figsize_from_map_extent(map_extent[region], base_height=6), layout='constrained')
                    ax = fig.add_subplot(1, 1, 1, projection=create_map_projection(region, central_longitudes=central_longitudes))
                    ax = setup_map_axes(ax, region, create_map_projection(region, central_longitudes=central_longitudes), map_extent=map_extent)

                    # Create scatter plot
                    sc = ax.scatter(lon_points, lat_points, c=data_points, 
                                    cmap=cmap, s=6, transform=ccrs.PlateCarree(),
                                    vmin=levels.min(), vmax=levels.max())
                    
                    if airstream != 'DI':
                        contour_50_map = ax.contour(lon, lat, percentage_grid, 
                                                    levels=[threshold_50], colors='red', alpha=1, linewidths=2,
                                                    transform=ccrs.PlateCarree(), linestyles='-', label='50%')
                        # contour_90_map = ax.contour(lon, lat, percentage_grid, 
                        #                             levels=[threshold_90], colors='red', linewidths=2,
                        #                             transform=ccrs.PlateCarree(), linestyles='-')
                        contour_95_map = ax.contour(lon, lat, percentage_grid, 
                                                    levels=[threshold_95], colors='darkred', linewidths=2,
                                                    transform=ccrs.PlateCarree(), linestyles='-', label='95%')
                    if airstream == 'DI':
                        if time_idx in range(0, 2):
                            contour_50_map = ax.contour(lon, lat, percentage_grid, 
                                                    levels=[threshold_50], colors='red', alpha=1, linewidths=2,
                                                    transform=ccrs.PlateCarree(), linestyles='-', label='50%')
                            # contour_90_map = ax.contour(lon, lat, percentage_grid, 
                            #                             levels=[threshold_90], colors='red', linewidths=2,
                            #                             transform=ccrs.PlateCarree(), linestyles='-')
                            contour_95_map = ax.contour(lon, lat, percentage_grid, 
                                                        levels=[threshold_95], colors='darkred', linewidths=2,
                                                        transform=ccrs.PlateCarree(), linestyles='-', label='95%')
                    # ax.clabel(contour_50_map, inline=True, fontsize=16, colors='red', fmt=lambda x: '50%')
                    # ax.clabel(contour_95_map, inline=True, fontsize=16, colors='darkred', fmt=lambda x: '95%')
                    # ax.set_extent(map_extent[region], crs=ccrs.PlateCarree())
                    # cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.05, aspect=50, extend='both')
                    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, aspect=50, extend='both', shrink=get_colorbar_shrink(map_extent[region], base_width=10.0, base_shrink=1.0))
                    cbar.set_label(units_dict[varname], size=20)
                    cbar.set_ticks(levels[::2])  # Use every second level for cleaner ticks
                    cbar.ax.tick_params(labelsize=16)
                    #ax.set_title(f"dTH (24hourly) at t={time_hours[time_idx]:.1f}h", fontsize=20)
                    legend_elements = [Line2D([0], [0], color='blue', lw=2, label='50% of trajectories'),
                                    Line2D([0], [0], color='red', lw=2, label='90% of trajectories'),
                                    Line2D([0], [0], color='darkred', lw=2, label='95% of trajectories')]
                    # ax.legend(handles=legend_elements, loc='upper right', fontsize=20)
                    if density_type == 'all':
                        plt.savefig(f"{fig_path}density_{airstream}_{region}_{timeframe}_{varname}_99.png", dpi=300, bbox_inches='tight')
                        plt.savefig(f"{fig_path}density_{airstream}_{region}_{timeframe}_{varname}_99.pdf", bbox_inches='tight')
                    else:
                        plt.savefig(f"{fig_path}density_{airstream}_{region}_{timeframe}_{varname}_{time_idx}.png", dpi=300, bbox_inches='tight')
                        plt.savefig(f"{fig_path}density_{airstream}_{region}_{timeframe}_{varname}_{time_idx}.pdf", bbox_inches='tight')
                    plt.show()

#%% Plot present and future of percentage grid
regions = ['JS', 'NA']
timeframes = ['present', 'future']
# airstreams = ['BCAR', 'DI', 'nonBCAR', 'LL', 'ALL']
airstreams = ['ALL', 'UL', 'LL']  # For now, only plot ALL airstream
airstreams = ['LL']  # For now, only plot ALL airstreams
# airstreams = ['ALL', 'BCAR', 'DI']  # For now, only plot BCAR airstream

#density files:
density_types = ['all']

# density_types = ['timesteps', 'all']
for density_type in density_types:
    if density_type == 'all':
        levels_dict['dTH_24'] = np.arange(0, 0.15, 0.025)  # Adjust levels for 24-hourly dTH
    else:
        levels_dict['dTH_24'] = np.arange(-3, 3.25, 0.25)
    for region in regions:
        for timeframe in timeframes:
            for airstream in airstreams:
                density_file = f"{path}density_{density_type}_CAO_{airstream}_trajs_{region}_{timeframe}_ALL_with_dTH"
                print(f"Processing density file: {density_file}")
                ds = load_density_data(density_file)
                lon = ds.longitude
                lat = ds.latitude

                if 'time' in ds.dims:
                    time_hours = ds.time.values / np.timedelta64(1, 'h')
                    print("Time coordinate in hours:", time_hours)
                
                # Create meshgrid for coordinates
                lon_2d, lat_2d = np.meshgrid(ds.longitude, ds.latitude)
                lon_flat = lon_2d.flatten()
                lat_flat = lat_2d.flatten()

                varname = 'dTH_24'
                cmap = cmc.oslo_r
                levels = levels_dict[varname]
                time_idxs = np.arange(len(ds.time)) if 'time' in ds.dims else [0]
                for time_idx in time_idxs:
                    #COUNTS
                    percentage_grid, threshold_50, threshold_90, threshold_95 = count_data(ds)

                    mask_95 = np.where(percentage_grid<=threshold_95, np.nan, 1)

                    # Apply mask to data
                    data = percentage_grid * mask_95

                    # Flatten the arrays for scatter plot
                    data_flat = data.flatten()

                    # Remove NaN values
                    mask = ~np.isnan(data_flat)
                    lon_points = lon_flat[mask]
                    lat_points = lat_flat[mask]
                    data_points = data_flat[mask]

                    if airstream in ['ALL', 'UL', 'LL']:
                        map_extent = map_extent
                        central_longitudes = central_longitudes
                    else:
                        map_extent = map_extent_DI
                        central_longitudes = central_longitudes_DI
                    fig = plt.figure(figsize=(10, 8), layout='constrained')
                    fig = plt.figure(figsize=get_figsize_from_map_extent(map_extent[region], base_height=6), layout='constrained')
                    ax = fig.add_subplot(1, 1, 1, projection=create_map_projection(region, central_longitudes=central_longitudes))
                    ax = setup_map_axes(ax, region, create_map_projection(region, central_longitudes=central_longitudes), map_extent=map_extent)

                    # Create scatter plot
                    sc = ax.scatter(lon_points, lat_points, c=data_points, 
                                    cmap=cmap, s=6, transform=ccrs.PlateCarree(), vmin = levels.min(), vmax = levels.max())
                    
                    if airstream != 'DI':
                        contour_50_map = ax.contour(lon, lat, percentage_grid, 
                                                    levels=[threshold_50], colors='red', alpha=1, linewidths=2,
                                                    transform=ccrs.PlateCarree(), linestyles='-', label='50%')
                        # contour_90_map = ax.contour(lon, lat, percentage_grid, 
                        #                             levels=[threshold_90], colors='red', linewidths=2,
                        #                             transform=ccrs.PlateCarree(), linestyles='-')
                        contour_95_map = ax.contour(lon, lat, percentage_grid, 
                                                    levels=[threshold_95], colors='darkred', linewidths=2,
                                                    transform=ccrs.PlateCarree(), linestyles='-', label='95%')
                    if airstream == 'DI':
                        if time_idx in range(0, 2):
                            contour_50_map = ax.contour(lon, lat, percentage_grid, 
                                                    levels=[threshold_50], colors='red', alpha=1, linewidths=2,
                                                    transform=ccrs.PlateCarree(), linestyles='-', label='50%')
                            # contour_90_map = ax.contour(lon, lat, percentage_grid, 
                            #                             levels=[threshold_90], colors='red', linewidths=2,
                            #                             transform=ccrs.PlateCarree(), linestyles='-')
                            contour_95_map = ax.contour(lon, lat, percentage_grid, 
                                                        levels=[threshold_95], colors='darkred', linewidths=2,
                                                        transform=ccrs.PlateCarree(), linestyles='-', label='95%')

                    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, aspect=50, extend='max', shrink=get_colorbar_shrink(map_extent[region], base_width=10.0, base_shrink=1.0))
                    cbar.set_label('[%]', size=20)
                    cbar.set_ticks(levels[::1])  # Use every level for cleaner ticks
                    cbar.ax.tick_params(labelsize=16)
                    legend_elements = [Line2D([0], [0], color='blue', lw=2, label='50% of trajectories'),
                                    Line2D([0], [0], color='red', lw=2, label='90% of trajectories'),
                                    Line2D([0], [0], color='darkred', lw=2, label='95% of trajectories')]
                    if density_type == 'all':
                        plt.savefig(f"{fig_path}density_percentage_{airstream}_{region}_{timeframe}_{varname}_99.png", dpi=300, bbox_inches='tight')
                        plt.savefig(f"{fig_path}density_percentage_{airstream}_{region}_{timeframe}_{varname}_99.pdf", bbox_inches='tight')
                    else:
                        plt.savefig(f"{fig_path}density_percentage_{airstream}_{region}_{timeframe}_{varname}_{time_idx}.png", dpi=300, bbox_inches='tight')
                        plt.savefig(f"{fig_path}density_percentage_{airstream}_{region}_{timeframe}_{varname}_{time_idx}.pdf", bbox_inches='tight')
                    plt.show()
#%%
# Plotting differences between present and future density data
regions = ['JS', 'NA']
airstreams = ['BCAR', 'DI', 'nonBCAR', 'LL', 'ALL']
airstreams = ['ALL', 'UL', 'LL']  # For now, only plot ALL airstream
airstreams = ['LL']
# airstreams = ['ALL']  # For now, only plot ALL airstreams

#density files:
# density_type = 'timesteps'
density_types = ['all']
# density_types = ['timesteps', 'all']

for density_type in density_types:
    for region in regions:
        for airstream in airstreams:
            if density_type == 'all':
                if airstream == 'ALL':
                    levels_dict['dTH_24'] = np.arange(-0.2, 0.22, 0.02)
                    levels_dict['dTH_24'] = np.arange(-0.5, 0.55, 0.05)
                else:
                    levels_dict['dTH_24'] = np.arange(-0.5, 0.55, 0.05)
            else:
                levels_dict['dTH_24'] = np.arange(-1, 1.1, 0.1)
            density_present_file = f"{path}density_{density_type}_CAO_{airstream}_trajs_{region}_present_ALL_with_dTH"
            density_future_file = f"{path}density_{density_type}_CAO_{airstream}_trajs_{region}_future_ALL_with_dTH"
            print(f"Processing present density file: {density_present_file}")
            ds_present = load_density_data(density_present_file)
            print(f"Processing future density file: {density_future_file}")
            ds_future = load_density_data(density_future_file)

            lon= ds_present.longitude
            lat = ds_present.latitude

            if 'time' in ds_present.dims:
                time_hours = ds_present.time.values / np.timedelta64(1, 'h')
                print("Time coordinate in hours:", time_hours)
            # Create meshgrid for coordinates
            lon_2d, lat_2d = np.meshgrid(ds_present.longitude, ds_present.latitude)
            lon_flat = lon_2d.flatten()
            lat_flat = lat_2d.flatten()

            varname = 'dTH_24'
            cmap = cmaps_dict[varname]
            levels = levels_dict[varname]
            time_idxs = np.arange(len(ds_present.time)) if 'time' in ds_present.dims else [0]

            for time_idx in time_idxs:
                # DATA
                data_present = mask_data(ds_present, varname, time_idx)
                data_future = mask_data(ds_future, varname, time_idx)

                count_present, threshold_50_present, threshold_90_present, threshold_95_present = count_data(ds_present)
                count_future, threshold_50_future, threshold_90_future, threshold_95_future = count_data(ds_future)
                
                mask_95_present = np.where(count_present <= threshold_95_present, np.nan, 1)
                mask_95_future = np.where(count_future <= threshold_95_future, np.nan, 1)

                data_present = data_present * mask_95_present
                data_future = data_future * mask_95_future

                data_diff = data_future - data_present

                # Flatten the arrays for scatter plot
                data_diff_flat = data_diff.flatten()

                # Remove NaN values
                mask_diff= ~np.isnan(data_diff_flat)
                lon_points_diff = lon_flat[mask_diff]
                lat_points_diff = lat_flat[mask_diff]
                data_points_diff = data_diff_flat[mask_diff]

                if airstream in ['ALL', 'BCAR', 'nonBCAR']:
                    map_extent = map_extent
                    central_longitudes = central_longitudes
                else:
                    map_extent = map_extent_DI
                    central_longitudes = central_longitudes_DI
                fig = plt.figure(figsize=(10, 8), layout='constrained')
                fig = plt.figure(figsize=get_figsize_from_map_extent(map_extent[region], base_height=6), layout='constrained')
                ax = fig.add_subplot(1, 1, 1, projection=create_map_projection(region, central_longitudes=central_longitudes))
                ax = setup_map_axes(ax, region, create_map_projection(region, central_longitudes=central_longitudes), map_extent=map_extent)

                # Create scatter plot for present
                sc_diff = ax.scatter(lon_points_diff, lat_points_diff, c=data_points_diff, 
                                        cmap=cmap, s=6, transform=ccrs.PlateCarree(),
                                        vmin=levels.min(), vmax=levels.max(),
                                        label='Difference')

                if airstream != 'DI':
                    contour_95_map_present = ax.contour(lon, lat, count_present, 
                                                        levels=[threshold_95_present], colors='darkred', linewidths=2,
                                                        transform=ccrs.PlateCarree(), linestyles='-')
                    contour_95_future = ax.contour(lon, lat, count_future, 
                                                        levels=[threshold_95_future], colors='darkred', linewidths=2,
                                                        transform=ccrs.PlateCarree(), linestyles='--')
                if airstream == 'DI':
                    if time_idx in range(0, 2):
                        contour_95_map_present = ax.contour(lon, lat, count_present, 
                                                            levels=[threshold_95_present], colors='darkred', linewidths=2,
                                                            transform=ccrs.PlateCarree(), linestyles='-')
                        contour_95_future = ax.contour(lon, lat, count_future, 
                                                        levels=[threshold_95_future], colors='darkred', linewidths=2,
                                                        transform=ccrs.PlateCarree(), linestyles='--')

                ax.set_extent(map_extent[region], crs=ccrs.PlateCarree())
                # cbar = plt.colorbar(sc_diff, ax=ax, orientation='vertical', pad=0.05, aspect=50, extend='both')
                cbar = plt.colorbar(sc_diff, ax=ax, orientation='horizontal', pad=0.05, aspect=50, extend='both', shrink=get_colorbar_shrink(map_extent[region], base_width=10.0, base_shrink=1.0))
                cbar.set_label('[Δ K 24h$^{-1}$]', size=20)
                cbar.set_ticks(levels[::2])
                cbar.ax.tick_params(labelsize=16)
                #ax.set_title(f"Difference in dTH (24hourly) at t={time_hours[time_idx]:.1f}h", fontsize=20)
                legend_elements = [Line2D([0], [0], color='blue', lw=2, label='50% of trajectories (present)'),
                                Line2D([0], [0], color='blue', lw=2, linestyle='--', label='50% of trajectories (future)')]
                # ax.legend(handles=legend_elements, loc='upper right', fontsize=20)
                if density_type == 'all':
                    plt.savefig(f"{fig_path}density_{airstream}_{region}_diff_{varname}_99.png", dpi=300, bbox_inches='tight')
                    plt.savefig(f"{fig_path}density_{airstream}_{region}_diff_{varname}_99.pdf", bbox_inches='tight')
                else:   
                    plt.savefig(f"{fig_path}density_{airstream}_{region}_diff_{varname}_{time_idx}.png", dpi=300, bbox_inches='tight')
                    plt.savefig(f"{fig_path}density_{airstream}_{region}_diff_{varname}_{time_idx}.pdf", bbox_inches='tight')
                plt.show()


#%%
density_file= f"{path}density_timesteps_CAO_BCAR_trajs_JS_present_ALL_with_dTH"

# Load density data
ds = load_density_data(density_file)
lon = ds.longitude
lat = ds.latitude

# Convert time coordinate from nanoseconds to hours
if 'time' in ds.coords:
    time_hours = ds.time.values / np.timedelta64(1, 'h')
    print("Time coordinate in hours:", time_hours)

# Create meshgrid for coordinates
lon_2d, lat_2d = np.meshgrid(ds.longitude, ds.latitude)
lon_flat = lon_2d.flatten()
lat_flat = lat_2d.flatten()

varname = 'dTH_24'
cmap = cmaps_dict[varname]
levels = levels_dict[varname]

time_idxs = np.arange(len(ds.time)) if 'time' in ds.dims else [0]



for time_idx in time_idxs:

    # DATA
    data = mask_data(ds, varname, time_idx)
    percentage_grid, threshold_50, threshold_90, threshold_95 = count_data(ds)

    # Flatten the arrays for scatter plot
    data_flat = data.flatten()

    # Remove NaN values
    mask = ~np.isnan(data_flat)
    lon_points = lon_flat[mask]
    lat_points = lat_flat[mask]
    data_points = data_flat[mask]

    fig = plt.figure(figsize=(15, 12), layout = 'constrained')
    ax = fig.add_subplot(1, 1, 1, projection=create_map_projection('JS'))
    ax = setup_map_axes(ax, 'JS', create_map_projection('JS'))

    # Create scatter plot
    sc = ax.scatter(lon_points, lat_points, c=data_points, 
                cmap=cmap, s=10, transform=ccrs.PlateCarree(),
                vmin=levels.min(), vmax=levels.max())
    
    contour_50_map = ax.contour(lon, lat, percentage_grid, 
                            levels=[threshold_50], colors='blue', linewidths=2,
                            transform=ccrs.PlateCarree(), linestyles='-')
    contour_90_map = ax.contour(lon, lat, percentage_grid, 
                                levels=[threshold_90], colors='red', linewidths=2,
                                transform=ccrs.PlateCarree(), linestyles='-')
    contour_95_map = ax.contour(lon, lat, percentage_grid, 
                                levels=[threshold_95], colors='darkred', linewidths=2,
                                transform=ccrs.PlateCarree(), linestyles='-')

    ax.set_extent(map_extent['JS'], crs=ccrs.PlateCarree())
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.05, aspect=50)
    cbar.set_label(units_dict[varname], size = 20)
    cbar.set_ticks(levels[::2])  # Use every other level for cleaner ticks
    cbar.ax.tick_params(labelsize=16)
    ax.set_title(f"dTH (6hourly) at t={time_hours[time_idx]:.1f}h", fontsize=20)
    legend_elements = [Line2D([0], [0], color='blue', lw=2, label='50% of trajectories'),
                    Line2D([0], [0], color='red', lw=2, label='90% of trajectories'),
                    Line2D([0], [0], color='darkred', lw=2, label='95% of trajectories')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=20)
    plt.savefig(f"{fig_path}density_DI_JS_present_{varname}_{time_hours[time_idx]:.0f}h.png", dpi=300, bbox_inches='tight')
    plt.show()

# COUNTS
counts = mask_data(ds, 'COUNT')
# Use count_data directly (it's already a numpy array)
count_data = counts

total_count = count_data

# Calculate percentage of total counts for each grid point
total_counts = np.nansum(total_count)  # Sum all counts, ignoring NaN values
percentage_grid = (total_count / total_counts) * 100

# Create a cumulative percentage grid by sorting grid points by their percentage
# and calculating cumulative sum
flat_percentages = percentage_grid.flatten()
flat_percentages = flat_percentages[~np.isnan(flat_percentages)]  # Remove NaN values
sorted_percentages = np.sort(flat_percentages)[::-1]  # Sort descending

# Calculate cumulative sum of percentages
cumsum_percentages = np.cumsum(sorted_percentages)

# Find threshold values for 50%, 90%, and 95% of total trajectories
threshold_50 = sorted_percentages[np.where(cumsum_percentages >= 50)[0][0]]
threshold_90 = sorted_percentages[np.where(cumsum_percentages >= 90)[0][0]]
threshold_95 = sorted_percentages[np.where(cumsum_percentages >= 95)[0][0]]


# Add contour lines for counts

contour_50_map = ax.contour(lon, lat, percentage_grid, 
                            levels=[threshold_50], colors='blue', linewidths=2,
                            transform=ccrs.PlateCarree(), linestyles='-')
contour_90_map = ax.contour(lon, lat, percentage_grid, 
                            levels=[threshold_90], colors='red', linewidths=2,
                            transform=ccrs.PlateCarree(), linestyles='-')
contour_95_map = ax.contour(lon, lat, percentage_grid, 
                            levels=[threshold_95], colors='darkred', linewidths=2,
                            transform=ccrs.PlateCarree(), linestyles='-')

ax.set_extent(map_extent['JS'], crs=ccrs.PlateCarree())
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.05, aspect=50)
cbar.set_label(units_dict[varname])
cbar.set_ticks(levels[::2])  # Use every other level for cleaner ticks
cbar.ax.tick_params(labelsize=10)
ax.set_title(f"Pressure at t={time_hours[0]:.1f}h", fontsize=14)
plt.show()

# %%
# Example usage
if __name__ == "__main__":
    
    # Plot overview for each region and all time steps
    for region in regions:
        density_file = f"{path}density_timesteps_CAO_BCAR_trajs_{region}_present_ALL_with_dTH"
        
        # Load data to check time dimension
        ds = load_density_data(density_file)
        if ds is None:
            continue
            
        if 'time' in ds.dims:
            # Plot for each time step
            for i in range(len(ds.time)):
                time_hours = timedelta_to_hours(ds.time.values[i])
                print(f"\nPlotting density overview for region: {region}, time: {time_hours:.1f}h")
                output_file = f"{fig_path}density_overview_{region}_t{time_hours:.0f}h.png"
                plot_density_overview(density_file, region=region, time_idx=i, output_file=output_file)
        else:
            # No time dimension, plot single overview
            print(f"\nPlotting density overview for region: {region}")
            output_file = f"{fig_path}density_overview_{region}.png"
            plot_density_overview(density_file, region=region, output_file=output_file)
    
    # Plot time evolution for specific variables
    for variable in ['COUNT', 'p', 'TH']:
        for region in regions:
            density_file = f"{path}density_timesteps_CAO_BCAR_trajs_{region}_present_ALL_with_dTH"
            print(f"\nPlotting time evolution of {variable} for region: {region}")
            output_file = f"{fig_path}time_evolution_{variable}_{region}.png"
            plot_time_evolution(density_file, variable, region=region, output_file=output_file)
    

# %%

"""Density maps of the 24-hourly potential temperature change along back-trajectories.

For the low-level air stream (LL) of each region, plots the 24-hourly change of the potential
temperature (dTH_24) on the Lagranto density grid, for the present and the future climate and
the difference between them. Grid points outside the area that holds 95% of the trajectories
are masked. The contours enclose 50% (red) and 95% (dark red) of the trajectories.

Usage (from any directory):

    python plot_density.py FIGURE [FIGURE ...]

where FIGURE is one of:
    F4_maps   present and future maps for JS and NA (Fig. 4a, b, d, e)
    F4_diff   future-present difference maps for JS and NA (Fig. 4c, f)

The density files are written by density_trajs.sh and read from $CAO_DENSITY_DIR. The
figures are written to $CAO_PLOTS_ROOT/density (see config/config.py).
"""

import argparse
import os
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
import config

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"

# Settings
path = f'{config.DENSITY_DIR}/'
fig_path = f'{config.PLOTS_ROOT}/density/'

airstream = 'LL'
varname = 'dTH_24'
unit_label = '[K 24h$^{-1}$]'
cmap = plt.get_cmap('RdBu_r')

# Color levels of the maps and of the difference maps
levels_maps = np.arange(-1, 1.1, 0.1)
levels_diff = np.arange(-0.5, 0.55, 0.05)

map_extent = {
    "JS": [60, 260, 25, 90],
    "NA": [130, 330, 25, 90]  # North Atlantic in 0-360 degrees east
}
central_longitudes = {
    "JS": 120,
    "NA": 270
}


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
    map_projection = ccrs.LambertAzimuthalEqualArea(central_longitude=central_longitudes[region], central_latitude=60)
    return map_projection


def setup_map_axes(ax, region, projection, map_extent):
    """Setup map axes with coastlines, borders, etc."""

    ax.add_feature(cfeature.LAND, facecolor=(217/256, 217/256, 217/256))
    ax.set_extent(map_extent[region], crs=ccrs.PlateCarree())
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


def get_colorbar_shrink(extent, base_width=12.0, base_shrink=0.6):
    """
    Calculates colorbar shrink based on a linear scaling formula.

    Parameters:
        extent: list
            Map extent [lon_min, lon_max, lat_min, lat_max]; the figure width follows from it
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


def count_data(ds, time_idx):
    """Percentage of all trajectories per grid point and the percentages that enclose
    50%, 90% and 95% of all trajectories."""
    total_count = mask_data(ds, 'COUNT', time_idx)

    # Calculate percentage of total counts for each grid point
    total_counts = np.nansum(total_count)  # Sum all counts, ignoring NaN values
    percentage_grid = (total_count / total_counts) * 100

    # Sort the grid points by their percentage and calculate the cumulative sum
    flat_percentages = percentage_grid.flatten()
    flat_percentages = flat_percentages[~np.isnan(flat_percentages)]  # Remove NaN values
    sorted_percentages = np.sort(flat_percentages)[::-1]  # Sort descending
    cumsum_percentages = np.cumsum(sorted_percentages)

    # Find threshold values for 50%, 90%, and 95% of total trajectories
    threshold_50 = sorted_percentages[np.where(cumsum_percentages >= 50)[0][0]]
    threshold_90 = sorted_percentages[np.where(cumsum_percentages >= 90)[0][0]]
    threshold_95 = sorted_percentages[np.where(cumsum_percentages >= 95)[0][0]]

    return percentage_grid, threshold_50, threshold_90, threshold_95


def density_file(region, timeframe):
    return f"{path}density_all_CAO_{airstream}_trajs_{region}_{timeframe}_ALL_with_dTH"


def new_map_axes(region):
    """Create a figure with a map of the region."""
    fig = plt.figure(figsize=get_figsize_from_map_extent(map_extent[region], base_height=6), layout='constrained')
    projection = create_map_projection(region, central_longitudes=central_longitudes)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax = setup_map_axes(ax, region, projection, map_extent=map_extent)
    return fig, ax


def add_colorbar(fig, ax, sc, region, label):
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, aspect=50, extend='both',
                        shrink=get_colorbar_shrink(map_extent[region], base_width=10.0, base_shrink=1.0))
    cbar.set_label(label, size=20)
    return cbar


def plot_maps(region, timeframe):
    """Map of dTH_24 for one region and timeframe (Fig. 4a, b, d, e)."""
    levels = levels_maps
    ds = load_density_data(density_file(region, timeframe))
    lon = ds.longitude
    lat = ds.latitude

    # Create meshgrid for coordinates
    lon_2d, lat_2d = np.meshgrid(ds.longitude, ds.latitude)
    lon_flat = lon_2d.flatten()
    lat_flat = lat_2d.flatten()

    time_idxs = np.arange(len(ds.time)) if 'time' in ds.dims else [0]
    for time_idx in time_idxs:
        data = mask_data(ds, varname, time_idx)
        percentage_grid, threshold_50, threshold_90, threshold_95 = count_data(ds, time_idx)

        # Mask the grid points outside the 95% contour of the trajectory counts
        mask_95 = np.where(percentage_grid <= threshold_95, np.nan, 1)
        data = data * mask_95

        # Flatten the arrays for scatter plot and remove NaN values
        data_flat = data.flatten()
        mask = ~np.isnan(data_flat)

        fig, ax = new_map_axes(region)
        sc = ax.scatter(lon_flat[mask], lat_flat[mask], c=data_flat[mask],
                        cmap=cmap, s=6, transform=ccrs.PlateCarree(),
                        vmin=levels.min(), vmax=levels.max())

        ax.contour(lon, lat, percentage_grid, levels=[threshold_50], colors='red', alpha=1, linewidths=2,
                   transform=ccrs.PlateCarree(), linestyles='-')
        ax.contour(lon, lat, percentage_grid, levels=[threshold_95], colors='darkred', linewidths=2,
                   transform=ccrs.PlateCarree(), linestyles='-')

        cbar = add_colorbar(fig, ax, sc, region, unit_label)
        cbar.set_ticks(levels[::2])  # Use every second level for cleaner ticks
        cbar.ax.tick_params(labelsize=16)

        plt.savefig(f"{fig_path}density_{airstream}_{region}_{timeframe}_{varname}_99.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{fig_path}density_{airstream}_{region}_{timeframe}_{varname}_99.pdf", bbox_inches='tight')
        plt.close(fig)


def plot_difference(region):
    """Future minus present map of dTH_24 for one region (Fig. 4c, f).
    The solid and dashed contours enclose 95% of the present and future trajectories."""
    levels = levels_diff
    ds_present = load_density_data(density_file(region, 'present'))
    ds_future = load_density_data(density_file(region, 'future'))

    lon = ds_present.longitude
    lat = ds_present.latitude
    lon_2d, lat_2d = np.meshgrid(ds_present.longitude, ds_present.latitude)
    lon_flat = lon_2d.flatten()
    lat_flat = lat_2d.flatten()

    time_idxs = np.arange(len(ds_present.time)) if 'time' in ds_present.dims else [0]
    for time_idx in time_idxs:
        data_present = mask_data(ds_present, varname, time_idx)
        data_future = mask_data(ds_future, varname, time_idx)

        count_present, _, _, threshold_95_present = count_data(ds_present, time_idx)
        count_future, _, _, threshold_95_future = count_data(ds_future, time_idx)

        # Mask the grid points outside the 95% contour of the trajectory counts
        data_present = data_present * np.where(count_present <= threshold_95_present, np.nan, 1)
        data_future = data_future * np.where(count_future <= threshold_95_future, np.nan, 1)
        data_diff_flat = (data_future - data_present).flatten()
        mask_diff = ~np.isnan(data_diff_flat)

        fig, ax = new_map_axes(region)
        sc_diff = ax.scatter(lon_flat[mask_diff], lat_flat[mask_diff], c=data_diff_flat[mask_diff],
                             cmap=cmap, s=6, transform=ccrs.PlateCarree(),
                             vmin=levels.min(), vmax=levels.max(),
                             label='Difference')

        ax.contour(lon, lat, count_present, levels=[threshold_95_present], colors='darkred', linewidths=2,
                   transform=ccrs.PlateCarree(), linestyles='-')
        ax.contour(lon, lat, count_future, levels=[threshold_95_future], colors='darkred', linewidths=2,
                   transform=ccrs.PlateCarree(), linestyles='--')

        ax.set_extent(map_extent[region], crs=ccrs.PlateCarree())
        cbar = add_colorbar(fig, ax, sc_diff, region, '[Δ K 24h$^{-1}$]')
        cbar.set_ticks(levels[::2])
        cbar.ax.tick_params(labelsize=16)

        plt.savefig(f"{fig_path}density_{airstream}_{region}_diff_{varname}_99.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{fig_path}density_{airstream}_{region}_diff_{varname}_99.pdf", bbox_inches='tight')
        plt.close(fig)


def f4_maps():
    for region in ('JS', 'NA'):
        for timeframe in ('present', 'future'):
            plot_maps(region, timeframe)


def f4_diff():
    for region in ('JS', 'NA'):
        plot_difference(region)


# Figure name -> function that produces it
FIGURES = {
    "F4_maps": f4_maps,
    "F4_diff": f4_diff,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("figures", nargs="+", choices=sorted(FIGURES), help="figures to produce")
    args = parser.parse_args()
    os.makedirs(fig_path, exist_ok=True)
    for name in args.figures:
        print(f"Producing {name}")
        FIGURES[name]()


if __name__ == "__main__":
    main()

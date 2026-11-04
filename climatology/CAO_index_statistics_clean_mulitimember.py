"""Ensemble statistics of marine cold-air outbreaks (CAOs) in the CESM2 large ensemble.

Produces the maps and bar charts of the CAO index (CAOI), the CAO frequency, the surface
evaporation and the contribution of the CAO index components, for the present and the
future climate, as used in the paper.

Usage (from any directory):

    python CAO_index_statistics_clean_mulitimember.py FIGURE [FIGURE ...]

where FIGURE is one of: F1, S1, S2, F2, S3 (see FIGURES below). The input data are read
from $CAO_DATA_ROOT and the figures are written to $CAO_PLOTS_ROOT/statistics/ensemble
(see config/config.py).
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm
import matplotlib
import pandas as pd
import numpy as np
from scipy.stats import norm
from matplotlib.patches import Patch
from scipy.stats import ttest_ind
import xarray as xr
import cftime
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import argparse
import gc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
import config

# Ensemble members, regions and timeframes are fixed by the study design.
member_ids = [1000, 1100, 1200, 1300, 1400]
fig_path = f'{config.PLOTS_ROOT}/statistics/ensemble'

# Plot style
plt.rcParams.update({
    "text.usetex": False
})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams.update({
    'font.size': 15.5,           # Default font size for ticks
    'axes.titlesize': 15.5,      # Title font size
    'axes.labelsize': 15.5,      # X and Y label font size
    'legend.fontsize': 15.5,     # Legend font size
    'xtick.labelsize': 10.34,     # X-axis tick font size
    'ytick.labelsize': 10.34      # Y-axis tick font size
})

# Map extents, box for the region mean and variable settings
map_extent = {
    "JS": [120, 160, 20.2617801, 64.554974],  # for map
    "NA": [270, 310, 20.2617801, 64.554974]   # for map
}

map_extent_precip = {
    "JS": [120, 160, 20.2617801, 64.554974],  # for map, use PlateCarree()
    "NA": [270, 310, 20.2617801, 64.554974]   # for map, use PlateCarree()
}

read_extent_precip = {
    "JS": [120-20, 260+20, 0, 90],  # for reading
    "NA": [230, 15+20, 0, 90]   # for reading old: 275-20

}

box = {
    "JS": [130, 150, 30, 60],  # for mask (50 north for Fig S3)
    "NA": [285, 295, 30, 45]   # for mask
}

VAR_CONFIG = {
    'CAOI': {
        'cmap_name': 'davos_r',
        'boundaries': np.arange(0, 19, 1),
        'cbarticks': np.arange(0, 19, 2),
        'labname': 'CAOI [K]'
    },
    'SHFLX_masked': {
        'cmap_name': 'viridis',
        'boundaries': np.arange(-200, 250, 50),
        'cbarticks': np.arange(-200, 300, 100),
        'labname': 'surface sensible heat flux [w m$^{-2}$]'
    },
    'LHFLX_masked': {
        'cmap_name': 'plasma',
        'boundaries': np.arange(-100, 150, 25),
        'cbarticks': np.arange(-100, 200, 50),
        'labname': 'surface latent heat flux [w m$^{-2}$]'
    },
    'SSTH_masked': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(273, 305, 2),
        'cbarticks': np.arange(273, 305, 4),
        'labname': '$\\Theta_{SST}$ [K]'
    },
    'TH_850_masked': {
        'cmap_name': 'plasma',
        'boundaries': np.arange(250, 290, 5),
        'cbarticks': np.arange(250, 300, 10),
        'labname': '$\\Theta_{850}$ [K]'
    },
    'PRECIP': {
        'cmap_name': 'BuPu',
        'boundaries': np.arange(0, 3.25, 0.25),
        'cbarticks': np.arange(0, 3.5, 0.5),
        'labname': 'precipitation [mm 6h$^{-1}$]'
    },
    'EVAP': {
        'cmap_name': 'YlGnBu',
        'boundaries': np.arange(0, 6.25, 0.5),
        'cbarticks': np.arange(0, 6.25, 1),
        'labname': 'evaporation [mm 6h$^{-1}$]'
    },
    'EP': {
        'cmap_name': 'RdYlBu_r',
        'boundaries': np.arange(-4, 4.5, 0.5),
        'cbarticks': np.arange(-4, 4.5, 1),
        'labname': 'net evaporation [mm 6h$^{-1}$]'
    },
    'SST': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(273, 305, 2),
        'cbarticks': np.arange(273, 305, 4),
        'labname': 'SST [K]'
    },
    'SSTH_ano': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(-10, 11, 1),
        'cbarticks': np.arange(-10, 11, 2),
        'labname': '$\\Theta_{SST}$ anomaly [K]'
    },
    'TH_850_ano': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(-15, 16, 1.5),
        'cbarticks': np.arange(-15, 16, 3),
        'labname': '$\\Theta_{850}$ anomaly [K]'
    }
}

VAR_CONFIG_DIFF = {
    'CAOI': {
        'cmap_name': 'vik',
        'boundaries': np.arange(-3, 3.5, 0.5),
        'cbarticks': np.arange(-3, 4, 1),
        'labname': 'Δ CAOI [K]'
    },
    'SHFLX_masked': {
        'cmap_name': 'RdBu_r',
        'boundaries': np.arange(-50, 55, 5),
        'cbarticks': np.arange(-50, 60, 10),
        'labname': 'Δ surface sensible heat flux [w m$^{-2}$]'
    },
    'LHFLX_masked': {
        'cmap_name': 'RdBu_r',
        'boundaries': np.arange(-30, 35, 5),
        'cbarticks': np.arange(-30, 40, 10),
        'labname': 'Δ surface latent heat flux [w m$^{-2}$]'
    },
    'SSTH_masked': {
        'cmap_name': 'RdBu_r',
        'boundaries': np.arange(-5, 6, 1),
        'cbarticks': np.arange(-5, 6, 1),
        'labname': 'Δ $\\Theta_{SST}$ [K]'
    },
    'TH_850_masked': {
        'cmap_name': 'RdBu_r',
        'boundaries': np.arange(-2, 2.5, 0.5),
        'cbarticks': np.arange(-2, 3, 1),
        'labname': 'Δ $\\Theta_{850}$ [K]'
    },
    'PRECIP': {
        'cmap_name': 'BrBG',
        'boundaries': np.arange(-0.8, 0.9, 0.1),
        'cbarticks': np.arange(-0.8, 0.9, 0.2),
        'labname': 'Δ precipitation [mm 6h$^{-1}$]'
    },
    'EVAP': {
        'cmap_name': 'vik',
        'boundaries': np.arange(-0.8, 0.9, 0.1),
        'cbarticks': np.arange(-0.8, 0.9, 0.2),
        'labname': 'Δ evaporation [mm 6h$^{-1}$]'
    },
    'EP': {
        'cmap_name': 'vik',
        'boundaries': np.arange(-1, 1.1, 0.2),
        'cbarticks': np.arange(-1, 1.1, 0.4),
        'labname': 'Δ net evaporation [mm 6h$^{-1}$]'
    },
    'SST': {
        'cmap_name': 'RdBu_r',
        'boundaries': np.arange(-5, 6, 1),
        'cbarticks': np.arange(-5, 6, 1),
        'labname': 'Δ SST [K]'
    },
    'TH_850_ano': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(-2, 2.2, 0.2),
        'cbarticks': np.arange(-2, 2.2, 0.4),
        'labname': 'Δ $\\Theta_{SST}$ anomaly [K]'
    },
    'SSTH_ano': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(-2, 2.2, 0.2),
        'cbarticks': np.arange(-2, 2.2, 0.4),
        'labname': '$\\Theta_{850}$ anomaly [K]'
    }

}

def djf_mask(ds, timeframe):
    """Create a mask for 10 DJF seasons in the dataset. Use as ds.sel(time=djf_mask(ds))"""
    # Assume ds is your original dataset
    time = ds['time'].to_index()
    if timeframe == "present":
        start_year_dec = 1990
        end_year_dec = 1999
        start_year_jan = start_year_dec + 1
        end_year_jan = end_year_dec +1
    elif timeframe == "future":
        start_year_dec = 2090
        end_year_dec = 2099
        start_year_jan = start_year_dec + 1
        end_year_jan = end_year_dec +1

    # Build a boolean mask for valid DJF times
    mask = (
        ((time.month == 12) & (start_year_dec <= time.year) & (time.year <= end_year_dec)) |
        ((np.isin(time.month, [1, 2])) & (start_year_jan <= time.year) & (time.year <= end_year_jan))
    )

    # Get the valid times
    valid_times = time[mask]
    return valid_times

def read_ensemble_data(region, timeframe, member_ids):
    ds_list = []
    for member_id in member_ids:
        path = f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/data_{timeframe}"
        ds = xr.open_dataset(path, decode_times=True, engine='netcdf4')
        ds = ds.expand_dims({'member': [member_id]})
        ds_list.append(ds)
    all_members = xr.concat(ds_list, dim='member')
    
    # Calculate event-box mean for this specific timeframe
    event_box_data = get_event_in_box_ensemble(all_members, region, timeframe, 'SSTH_masked', member_ids)
    weights = np.cos(np.deg2rad(event_box_data.lat))
    weights.name = "weights"
    
    ssth_event_box_mean = event_box_data.weighted(weights).mean(dim=['time', 'lat', 'lon'], skipna=True)
    
    event_box_data_th = get_event_in_box_ensemble(all_members, region, timeframe, 'TH_850_masked', member_ids)
    th850_event_box_mean = event_box_data_th.weighted(weights).mean(dim=['time', 'lat', 'lon'], skipna=True)
    
    # Calculate anomalies
    ssth_ano = all_members.SSTH_masked - ssth_event_box_mean
    ssth_ano.attrs["units"] = "K"
    all_members["SSTH_ano"] = ssth_ano
    
    th850_ano = all_members.TH_850_masked - th850_event_box_mean
    th850_ano.attrs["units"] = "K"
    all_members["TH_850_ano"] = th850_ano
    
    return all_members

def read_ensemble_data_precip(region, timeframe, member_ids):
    """
    Read ensemble data for precipitation, evaporation, and P-E.
    Returns: ds_precip with dims (time, lat, lon, member)
    """
    ds_list = []
    for member_id in member_ids:
        path = f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}"
        ds_precip = xr.open_dataset(path, decode_times=True, engine='netcdf4')
        ds_precip = ds_precip.expand_dims({'member': [member_id]})
        # Handle longitude selection for regions that cross 0° meridian
        if read_extent_precip[region][1] < read_extent_precip[region][0]:
            # For NA: 275° to 30° crosses 0° meridian
            # Select 275° to 360° and 0° to 30°, then concatenate
            ds_precip_west = ds_precip.sel(
                lon=slice(read_extent_precip[region][0], 360), 
                lat=slice(read_extent_precip[region][2],read_extent_precip[region][3])
            )
            ds_precip_east = ds_precip.sel(
                lon=slice(0, read_extent_precip[region][1]), 
                lat=slice(read_extent_precip[region][2], read_extent_precip[region][3])
            )
            ds_precip = xr.concat([ds_precip_west, ds_precip_east], dim='lon')
        else:
            # Normal case
            ds_precip = ds_precip.sel(
                lon=slice(read_extent_precip[region][0], read_extent_precip[region][1]), 
                lat=slice(read_extent_precip[region][2], read_extent_precip[region][3])
                    # map_extent_precip[region][2],
                        #   map_extent_precip[region][3]
            )        
        ds_list.append(ds_precip)
    all_members= xr.concat(ds_list, dim='member')
    # Calculate PRECIP as mm/6h and add to dataset
    precip = all_members.PRECT_masked * (1000 * 21600)
    precip.attrs["units"] = "mm/6h"
    all_members["PRECIP"] = precip
    #Calculate EVAP as mm/6h and add to dataset
    evap = (all_members.LHFLX_masked * 1000 * 21600)/(1000 * 2.45e6)  # Convert W/m² to mm/6h
    evap.attrs["units"] = "mm/6h"
    all_members["EVAP"] = evap
    #Calculate P-E as mm/6h and add to dataset
    pe = all_members.EVAP - all_members.PRECIP
    pe.attrs["units"] = "mm/6h"
    all_members["EP"] = pe
    return all_members

def read_ensemble_data_ice(region, timeframe, member_ids):
    ds_list = []
    for member_id in member_ids:
        path = f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/data_ice_{timeframe}"
        ds = xr.open_dataset(path, decode_times=True, engine='netcdf4')
        ds = ds.expand_dims({'member': [member_id]})
        ds_list.append(ds)
    all_members= xr.concat(ds_list, dim='member')
    return all_members

def calculate_space_mean(CAOI_djf):
    """(time,lat, lon) --> (time,) with latitudinal weights"""
    # use latitudinal weights to calculate the mean
    weights = np.cos(np.deg2rad(CAOI_djf.lat))
    weights.name = "weights"
    CAOI_djf_weighted = CAOI_djf.weighted(weights)
    CAOI_djf__mean_space = CAOI_djf_weighted.mean(dim=["lat", "lon"], skipna=True)

    #method from evolution script
    # CAOI_djf__mean_space = np.nanmean(np.nanmean(CAOI_djf, axis=1), axis=1)

    return CAOI_djf__mean_space

def calculate_space_mean_no_weights(CAOI_djf):
    """(time,lat, lon) --> (time,) without latitudinal weights"""
    CAOI_djf__mean_space = CAOI_djf.mean(dim=["lat", "lon"], skipna=True)
    return CAOI_djf__mean_space

def calculate_time_mean(CAOI_djf):
    """(time,lat, lon) --> (lat, lon)"""
    # use latitudinal weights to calculate the mean
    weights = np.cos(np.deg2rad(CAOI_djf.lat))
    weights.name = "weights"
    CAOI_djf_weighted = CAOI_djf.weighted(weights)
    CAOI_djf__mean_time = CAOI_djf.mean(dim = ['time'], skipna=True, keep_attrs=True)
    return CAOI_djf__mean_time

def event_mask_time(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    """
    file_path = f"{config.DATA_ROOT}/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/data_{timeframe}", decode_times=True)
    time_coord = ds.time

    # Convert summary_events["time"] to the same type as ds.time
    if np.issubdtype(time_coord.dtype, np.datetime64):
        summary_events["time"] = pd.to_datetime(summary_events["time"], format="%Y-%m-%d-%H")
    else:
        import cftime
        summary_events["time"] = summary_events["time"].apply(
            lambda x: cftime.DatetimeNoLeap(*map(int, x.split("-")))
        )

    mask = time_coord.isin(summary_events["time"])
    return mask

def event_plus0_mask_time(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    This mask includes the event time only (no additional days).
    """
    file_path = f"{config.DATA_ROOT}/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
    time_coord = ds.time

    # Convert summary_events["time"] to the same type as ds.time
    if np.issubdtype(time_coord.dtype, np.datetime64):
        summary_events["time"] = pd.to_datetime(summary_events["time"], format="%Y-%m-%d-%H")
    else:
        import cftime
        summary_events["time"] = summary_events["time"].apply(
            lambda x: cftime.DatetimeNoLeap(*map(int, x.split("-")))
        )

    # Create a mask that includes the event time and +6h, +12h, -6h, -12h
    # This is done by checking if the time coordinate is within the event time ± 6 hours
    time_plus_6h = summary_events["time"] + pd.Timedelta(hours= 6)
    time_minus_6h = summary_events["time"] - pd.Timedelta(hours= 6)
    time_plus_12h = summary_events["time"] + pd.Timedelta(hours= 12)
    # time_minus_12h = summary_events["time"] - pd.Timedelta(hours= 12)
    mask = (
        time_coord.isin(summary_events["time"]) |
        time_coord.isin(time_plus_6h) |
        time_coord.isin(time_minus_6h) |
        time_coord.isin(time_plus_12h)
    )

    return mask

def event_plus1_mask_time(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    This mask includes the event time plus one day.
    """
    file_path = f"{config.DATA_ROOT}/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
    time_coord = ds.time

    # Convert summary_events["time"] to the same type as ds.time
    if np.issubdtype(time_coord.dtype, np.datetime64):
        summary_events["time"] = pd.to_datetime(summary_events["time"], format="%Y-%m-%d-%H")
    else:
        import cftime
        summary_events["time"] = summary_events["time"].apply(
            lambda x: cftime.DatetimeNoLeap(*map(int, x.split("-")))
        )

    # Create a mask that includes the event time plus one day
    time_0 = summary_events["time"] + pd.Timedelta(days=1)
    time_plus_6h = summary_events["time"] + pd.Timedelta(hours=6) + pd.Timedelta(days=1)
    time_minus_6h = summary_events["time"] - pd.Timedelta(hours=6) + pd.Timedelta(days=1)
    time_plus_12h = summary_events["time"] + pd.Timedelta(hours=12) + pd.Timedelta(days=1)
    # time_minus_12h = summary_events["time"] - pd.Timedelta(hours=12) + pd.Timedelta(days=1)
    mask = (
        time_coord.isin(time_0) |
        time_coord.isin(time_plus_6h) |
        time_coord.isin(time_minus_6h) |
        time_coord.isin(time_plus_12h)
    )
    return mask

def event_plus2_mask_time(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    This mask includes the event time plus two days.
    """
    file_path = f"{config.DATA_ROOT}/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"{config.DATA_ROOT}/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
    time_coord = ds.time

    # Convert summary_events["time"] to the same type as ds.time
    if np.issubdtype(time_coord.dtype, np.datetime64):
        summary_events["time"] = pd.to_datetime(summary_events["time"], format="%Y-%m-%d-%H")
    else:
        import cftime
        summary_events["time"] = summary_events["time"].apply(
            lambda x: cftime.DatetimeNoLeap(*map(int, x.split("-")))
        )

    # Create a mask that includes the event time plus two days
    time_0 = summary_events["time"] + pd.Timedelta(days=2)
    time_plus_6h = summary_events["time"] + pd.Timedelta(hours=6) + pd.Timedelta(days=2)
    time_minus_6h = summary_events["time"] - pd.Timedelta(hours=6) + pd.Timedelta(days=2)
    time_plus_12h = summary_events["time"] + pd.Timedelta(hours=12) + pd.Timedelta(days=2)
    # time_minus_12h = summary_events["time"] - pd.Timedelta(hours=12) + pd.Timedelta(days=2)
    mask = (
        time_coord.isin(time_0) |
        time_coord.isin(time_plus_6h) |
        time_coord.isin(time_minus_6h) |
        time_coord.isin(time_plus_12h)
    )
    return mask

def get_masked_ensemble(ds_ens, region, timeframe, varname, mask_func, member_ids):
    """ Use for event_mask_time or percentile_mask_time, returns DataArray with dims (time, lat, lon)"""
    masked_list = []
    for i, member_id in enumerate(member_ids):
        ds_member = ds_ens.sel(member=member_id)
        mask = mask_func(region, timeframe, member_id)
        masked = ds_member[varname].sel(time=mask)
        masked_list.append(masked)
    
    # Align and stack along 'member'
    # aligned_list = xr.align(*masked_list, join = 'outer')
    # return xr.concat(aligned_list, dim='member')
    member_data = []
    for i, da in enumerate(masked_list):  # assuming each item in datasets is a DataArray
        n_time = da.sizes['time']
        
        # Save original time values
        time_vals = da['time'].values

        # Drop time as coordinate by resetting its role as a dimension
        da = da.copy()
        da = da.drop_vars('time')  # remove time as coordinate
        da.coords['time_value'] = ('time', time_vals)  # re-add time as variable

        # Add member ID
        da.coords['member'] = ('time', [i] * n_time)

        member_data.append(da)

    # Concatenate along time
    return xr.concat(member_data, dim='time')  # no time merging now!

def get_event_in_box_ensemble(ds_ens, region, timeframe, varname, member_ids):
    """Returns DataArray with dims (time, lat, lon) for all ensemble members in the box."""
    ds_box_masked = ds_ens.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3]))
    ds_time_masked = get_masked_ensemble(ds_box_masked, region, timeframe, varname, event_mask_time, member_ids)
    return ds_time_masked

def get_djf_in_box_ensemble(ds_ens, region, timeframe, varname, member_ids):
    """Returns DataArray with dims (time, lat, lon) for all DJF times and all ensemble members in the box."""
    # Select box
    ds_box = ds_ens.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3]))
    # Select only DJF times for each member
    masked_list = []
    for member_id in member_ids:
        ds_member = ds_box.sel(member=member_id)
        # DJF mask for this member
        djf_times = djf_mask(ds_member, timeframe)
        mask = ds_member.time.isin(djf_times)
        masked = ds_member[varname].sel(time=mask)
        masked_list.append(masked)
    # Concatenate along time
    member_data = []
    for i, da in enumerate(masked_list):
        n_time = da.sizes['time']
        time_vals = da['time'].values
        da = da.copy()
        da = da.drop_vars('time')
        da.coords['time_value'] = ('time', time_vals)
        da.coords['member'] = ('time', [i] * n_time)
        member_data.append(da)
    return xr.concat(member_data, dim='time')

def cmap_pcolormesh(cmap_name, boundaries, one_way=True):
    # Check if it's a cmcrameri colormap first
    if hasattr(cmc, cmap_name):
        base_cmap = getattr(cmc, cmap_name)
        # cmcrameri changed the direction of "bilbao" between versions 1.2 and 1.7. The figures use the
        # version that runs from light (low values) to dark red (high values).
        if cmap_name == 'bilbao' and sum(base_cmap(0.0)[:3]) < sum(base_cmap(1.0)[:3]):
            base_cmap = base_cmap.reversed()
    else:
        # Fall back to matplotlib colormaps
        base_cmap = plt.cm.get_cmap(cmap_name, len(boundaries) + 1)
    
    if hasattr(base_cmap, 'resampled'):
        discrete_cmap = base_cmap.resampled(len(boundaries) + 1)
    else:
        colors = base_cmap(np.linspace(0, 1, len(boundaries) + 1))
        discrete_cmap = ListedColormap(colors)
    
    colors = list(discrete_cmap(np.linspace(0, 1, len(boundaries) + 1)))
    
    under_color = colors[0]
    over_color = colors[-1]
    main_colors = colors[1:-1]

    cmap = ListedColormap(main_colors)
    cmap.set_under(under_color)
    cmap.set_over(over_color)

    return cmap, len(main_colors)

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

def plot_single_map(data_2d, lons, lats, region, timeframe_name, extent, cmap_name, boundaries, varname, labname = None, cbarticks = None, extend='neither', save_path=False, lag = 0, significance_mask = None):
    # set up map projection, event mask and extent depending on variable and region
    if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
        # For precipitation, evaporation, and P-E, use the map extent for the box
        if region == "NA":
            # For North Atlantic, use a different projection
            map_projection = ccrs.AzimuthalEquidistant(central_longitude=-55, central_latitude=50)
            map_projection = ccrs.PlateCarree()
            # extent = [map_extent_precip[region][0]-360, map_extent_precip[region][1]-15, map_extent_precip[region][2], map_extent_precip[region][3]]
            extent = map_extent_precip[region]
        elif region == 'JS':
            map_projection = ccrs.AzimuthalEquidistant(central_longitude=160, central_latitude=50)
            map_projection = ccrs.PlateCarree()
            extent = map_extent_precip[region]
        
        if lag == 0:
            event_mask_time_func = event_plus0_mask_time
        elif lag == 1:
            event_mask_time_func = event_plus1_mask_time
        elif lag == 2:
            event_mask_time_func = event_plus2_mask_time
        area_slice = read_extent_precip[region] # to match lat/lon of sea ice data with precipitation data

    else:
        # For other variables, use the box defined in the VAR_CONFIG
        map_projection = ccrs.PlateCarree()
        extent = map_extent[region]
        event_mask_time_func = event_mask_time
        area_slice = map_extent[region] # to match sea ice data with CAOI data

    if timeframe_name == "10th_present" or timeframe_name == 'present':
        ds_name = 'present'
        # SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SSTH_masked", event_mask_time, member_ids)
            time_slice = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
        else:
            time_slice = djf_mask(data_dict_ice[f"{region}_present"], 'present')
        
        print(f"Using time slice for {region} {ds_name}: {time_slice}")

        ice_time_sliced = data_dict_ice[f"{region}_present"].sel(
            time = time_slice)
        ice_masked_ensemble = get_masked_ensemble(ice_time_sliced, region, 'present', 'ICEFRAC', event_mask_time_func, member_ids)
        
        print('area_slice:', area_slice)
        if area_slice[1] < area_slice[0]:
            # For NA: 275° to 30° crosses 0° meridian
            # Select 275° to 360° and 0° to 30°, then concatenate
            sea_ice_west = ice_masked_ensemble.sel(
                lon=slice(area_slice[0], 360), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice_east = ice_masked_ensemble.sel(
                lon=slice(0, area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice = xr.concat([sea_ice_west, sea_ice_east], dim='lon')
            calculate_time_mean(sea_ice_west).plot()
            calculate_time_mean(sea_ice_east).plot()
        else:
            # Normal case
            sea_ice = ice_masked_ensemble.sel(
                lon=slice(area_slice[0], area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
            calculate_time_mean(sea_ice).plot() 
        print('Sea ice data shape:', sea_ice.shape)
        print(sea_ice.lat.max(), sea_ice.lat.min(), sea_ice.lon.max(), sea_ice.lon.min())
        

    elif timeframe_name == "10th_future" or timeframe_name == 'future':
        ds_name = 'future'
        # SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SSTH_masked", event_mask_time, member_ids)
            time_slice = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
        else:
            time_slice = djf_mask(data_dict_ice[f"{region}_future"], ds_name)
        
        print(f"Using time slice for {region} {ds_name}: {time_slice}")

        ice_time_sliced = data_dict_ice[f"{region}_future"].sel(
            time = time_slice)
        ice_masked_ensemble = get_masked_ensemble(ice_time_sliced, region, 'future', 'ICEFRAC', event_mask_time_func, member_ids)
        
        if area_slice[1] < area_slice[0]:
            # For NA: 275° to 30° crosses 0° meridian
            # Select 275° to 360° and 0° to 30°, then concatenate
            sea_ice_west = ice_masked_ensemble.sel(
                lon=slice(area_slice[0], 360), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice_east = ice_masked_ensemble.sel(
                lon=slice(0, area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice = xr.concat([sea_ice_west, sea_ice_east], dim='lon')
            calculate_time_mean(sea_ice_west).plot()
            calculate_time_mean(sea_ice_east).plot()
        else:
            # Normal case
            sea_ice = ice_masked_ensemble.sel(
                lon=slice(area_slice[0], area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
            calculate_time_mean(sea_ice).plot()
 
        print('Sea ice data shape:', sea_ice.shape)
        print(sea_ice.lat.max(), sea_ice.lat.min(), sea_ice.lon.max(), sea_ice.lon.min())

    elif timeframe_name == "10th_diff" or timeframe_name == 'diff':
        # SST_present = get_masked_ensemble(data_dict[f"{region}_present"], region, "present", "SST", event_mask_time_func, member_ids)
        # SST_future = get_masked_ensemble(data_dict[f"{region}_future"], region, "future", "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SST_present = get_masked_ensemble(data_dict[f"{region}_present"], region, "present", "SSTH_masked", event_mask_time, member_ids)
            SST_present = SST_present.where(SST_present > 273)  # Mask out SST values below 0°C
            SST_future = get_masked_ensemble(data_dict[f"{region}_future"], region, "future", "SSTH_masked", event_mask_time, member_ids)
            SST_future = SST_future.where(SST_future > 273)  # Mask out SST values below 0°C
            time_slice_present = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
            time_slice_future = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
        else:
            time_slice_present = djf_mask(data_dict_ice[f"{region}_present"], "present")
            time_slice_future = djf_mask(data_dict_ice[f"{region}_future"], "future")

        ice_present_time_sliced = data_dict_ice[f"{region}_present"].sel(
            time = time_slice_present)
        ice_present_masked_ensemble = get_masked_ensemble(ice_present_time_sliced, region, 'present', 'ICEFRAC', event_mask_time_func, member_ids)
        if area_slice[1] < area_slice[0]:
            # For NA: 275° to 30° crosses 0° meridian
            # Select 275° to 360° and 0° to 30°, then concatenate
            sea_ice_present_west = ice_present_masked_ensemble.sel(
                lon=slice(area_slice[0], 360), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice_present_east = ice_present_masked_ensemble.sel(
                lon=slice(0, area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice_present = xr.concat([sea_ice_present_west, sea_ice_present_east], dim='lon')
        else:
            # Normal case
            sea_ice_present = ice_present_masked_ensemble.sel(
                lon=slice(area_slice[0], area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
        print('Sea ice present data shape:', sea_ice_present.shape)
        print(sea_ice_present.lat.max(), sea_ice_present.lat.min(), sea_ice_present.lon.max(), sea_ice_present.lon.min())
        
        ice_future_time_sliced = data_dict_ice[f"{region}_future"].sel(
            time = time_slice_future)
        ice_future_masked_ensemble = get_masked_ensemble(ice_future_time_sliced, region, 'future', 'ICEFRAC', event_mask_time_func, member_ids)
        if area_slice[1] < area_slice[0]:
            # For NA: 275° to 30° crosses 0° meridian
            # Select 275° to 360° and 0° to 30°, then concatenate
            sea_ice_future_west = ice_future_masked_ensemble.sel(
                lon=slice(area_slice[0], 360), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice_future_east = ice_future_masked_ensemble.sel(
                lon=slice(0, area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
            sea_ice_future = xr.concat([sea_ice_future_west, sea_ice_future_east], dim='lon')
        else:
            # Normal case
            sea_ice_future = ice_future_masked_ensemble.sel(
                lon=slice(area_slice[0], area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
        print('Sea ice future data shape:', sea_ice_future.shape)
        print(sea_ice_future.lat.max(), sea_ice_future.lat.min(), sea_ice_future.lon.max(), sea_ice_future.lon.min())

    # fig, ax = plt.subplots(1, 1, figsize=(15, 12), subplot_kw=dict(projection=map_projection), layout='constrained')
    #fontsize=20, 16, shrink = 0.6
    # fig, ax = plt.subplots(1, 1, figsize=(10, 8), subplot_kw=dict(projection=map_projection), layout='constrained')
    #fontsize = 25, 20, shrink = 0.6
    fig, ax = plt.subplots(1, 1, figsize=get_figsize_from_map_extent(extent, base_height=7), subplot_kw=dict(projection=map_projection), layout='constrained')

    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.coastlines(color='black')
    ax.set_facecolor("white")

    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 15, "color": "k"}
    gl.ylabel_style = {"size": 15, "color": "k"}
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(0, 100, 10))

    cmap, n_main = cmap_pcolormesh(cmap_name, boundaries) #main colors, without extended colors
    print('cmap.colors', np.shape(cmap.colors))
    print('n_main:', n_main)

    if timeframe_name == "10th_diff" or timeframe_name == 'diff':

        # Convert main colors to array and set alpha
        main_colors = np.array(cmap.colors)
        main_colors[:, -1] = 0.8  # Set alpha to 0.8 for main colors
        
        # Create new colormap with modified alpha
        cmap_alpha = ListedColormap(main_colors)
        
        # Handle extension colors with alpha
        if extend in ['max', 'both']:
            over_color = list(cmap.get_over())
            if len(over_color) == 3:  # RGB
                over_color.append(0.8)  # Add alpha
            elif len(over_color) == 4:  # RGBA
                over_color[3] = 0.8  # Modify alpha
            cmap_alpha.set_over(over_color)
        
        if extend in ['min', 'both']:
            under_color = list(cmap.get_under())
            if len(under_color) == 3:  # RGB
                under_color.append(0.8)  # Add alpha
            elif len(under_color) == 4:  # RGBA
                under_color[3] = 0.8  # Modify alpha
            cmap_alpha.set_under(under_color)
        
        cmap = cmap_alpha




    ##### test
    #####
    norm = BoundaryNorm(boundaries, ncolors=n_main, clip=False)
    if varname in ['EP', 'SSTH_ano', 'TH_850_ano', 'CAOI']:
        extend = 'both'  # P-E can be negative and positive
    
    if varname in [
        # 'PRECIP', 'EVAP', 
        'SST']:
        # Convert longitudes to [-180, 180] range for plotting
        if region == "NA":
            # Convert longitudes > 180 to negative values
            lons_plot = lons.where(lons <= 180, lons - 360)
            # Sort by longitude to ensure monotonic coordinates
            sort_idx = np.argsort(lons_plot)
            lons_plot = lons_plot[sort_idx]
            data_2d = data_2d[:, sort_idx]
            # Also sort significance mask if provided
            if significance_mask is not None:
                significance_mask = significance_mask[:, sort_idx]
        else:
            lons_plot = lons
        im = ax.contourf(lons_plot, lats, data_2d, levels=boundaries, 
                     transform=ccrs.PlateCarree(), cmap=cmap, extend=extend)
    else:
        lons_plot = lons
        im = ax.pcolormesh(lons_plot, lats, data_2d, shading='auto',
                       transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)

    lons_SST = data_dict[f"{region}_present"].lon
    lats_SST = data_dict[f"{region}_present"].lat
    
    if timeframe_name == '10th_present' or timeframe_name == '10th_future' or timeframe_name == 'present' or timeframe_name == 'future':
        # ax.contour(lons_SST, lats_SST, calculate_time_mean(SST), levels=[273], colors='green', linewidths=3, transform=ccrs.PlateCarree(), alpha=1)
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SST_mean = calculate_time_mean(SST)
            # Plot contour lines (e.g., every 2K from 273 to 305)
            dSST_dlat = SST_mean.differentiate('lat')
            dSST_dlon = SST_mean.differentiate('lon')

            # Calculate gradient magnitude
            grad_mag = np.sqrt(dSST_dlat**2 + dSST_dlon**2)

            # Create mask where gradient magnitude < 1
            mask = grad_mag < 100

            # Mask SST values
            print('Mean over box (graphically): ', calculate_space_mean(SST_mean.sel(lat=slice(box[region][2], box[region][3]), lon=slice(box[region][0], box[region][1]))))
            event_SSTH_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
            event_SSTH_present_space_mean = calculate_space_mean(event_SSTH_present)
            print('Mean over box for event SSTH_present: ', event_SSTH_present_space_mean.values.flatten().mean())
            SST_mean = SST_mean.where(mask).sel(lat=slice(extent[2], extent[3]), lon=slice(extent[0], extent[1]))
            
            SST_contour = ax.contour(SST_mean.lon, SST_mean.lat, SST_mean, levels = np.arange(275, 305, 5), colors = 'black', linewidths = 2, transform=ccrs.PlateCarree())
            ax.clabel(SST_contour, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
            # ax.contour(lons_plot, lats, calculate_time_mean(SLP), levels=np.arange(-6, 0, 2), colors='dimgrey', linestyles = 'dashed',linewidths=2, transform=ccrs.PlateCarree(), alpha=1)
            # ax.contour(lons_plot, lats, calculate_time_mean(SLP), levels=np.arange(0, 8, 2), colors='dimgrey', linewidths=2, transform=ccrs.PlateCarree(), alpha=1)
        if varname == 'CAOI':
            # Get SST masked with the same event mask and extent as CAOI
            SST_event_masked = get_masked_ensemble(
                data_dict[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time, member_ids
            )
            # Select the same extent as CAOI
            SST_event_masked = SST_event_masked.sel(
                lon=slice(extent[0], extent[1]), lat=slice(extent[2], extent[3])
            )
            # Calculate the time mean for contour plotting
            SST_mean = calculate_time_mean(SST_event_masked)
            # Plot contour lines (e.g., every 2K from 273 to 305)
            dSST_dlat = SST_mean.differentiate('lat')
            dSST_dlon = SST_mean.differentiate('lon')

            # Calculate gradient magnitude
            grad_mag = np.sqrt(dSST_dlat**2 + dSST_dlon**2)

            # Create mask where gradient magnitude < 1
            mask = grad_mag < 100

            # Mask SST values
            SST_mean = SST_mean.where(mask)
            # SST_contour = ax.contour(
            #     SST_mean.lon, SST_mean.lat, SST_mean,
            #     levels=np.arange(275, 310, 5),
            #     colors='black', linewidths=2,
            #     transform=ccrs.PlateCarree()
            # )
            # ax.clabel(SST_contour, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
        
        print("lons_plot range:", lons_plot.min().values, "to", lons_plot.max().values)
        print("sea_ice.lon range:", sea_ice.lon.min().values, "to", sea_ice.lon.max().values)
        print("Coordinates match:", np.array_equal(lons_plot.values, sea_ice.lon.values))
        # Check if there's actually sea ice above the threshold before plotting
        sea_ice_mean = calculate_time_mean(sea_ice)
        max_ice_concentration = sea_ice_mean.max().values
        
        print(f"Max sea ice concentration for {timeframe_name}: {max_ice_concentration}")
        
        if max_ice_concentration >= 0.5:
            ax.contour(lons_plot, sea_ice_mean.lat, sea_ice_mean,
                    levels=[0.5], colors='green', linewidths=2,
                    transform=ccrs.PlateCarree(), alpha=1)
            print(f"Drew 50% sea ice contour for {timeframe_name}")
        else:
            print(f"No 50% sea ice contour available for {timeframe_name} (max: {max_ice_concentration:.3f})")
    elif timeframe_name == '10th_diff' or timeframe_name == 'diff':
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SST_diff_mean = calculate_time_mean(SST_future)-calculate_time_mean(SST_present)
            SST_diff_mean = SST_diff_mean.sel(lat=slice(extent[2], extent[3]), lon=slice(extent[0], extent[1]))
            SST_contour_neg = ax.contour(SST_diff_mean.lon, SST_diff_mean.lat, SST_diff_mean, levels = [-4, -3, -2, -1], colors = 'black', negative_linestyles = 'dashed', linewidths = 2, transform=ccrs.PlateCarree())
            SST_contour_pos = ax.contour(SST_diff_mean.lon, SST_diff_mean.lat, SST_diff_mean, levels = [1, 2, 3, 4, 5, 6], colors = 'black', negative_linestyles = 'dashed', linewidths = 2, transform=ccrs.PlateCarree())
            ax.clabel(SST_contour_neg, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
            ax.clabel(SST_contour_pos, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
        # if varname == 'CAOI':
        #     # Get SST masked with the same event mask and extent as CAOI present and future
        #     SST_present_event_masked = get_masked_ensemble(
        #         data_dict[f"{region}_present"], region, "present", "SST", event_mask_time, member_ids
        #     )
        #     SST_present_event_masked = SST_present_event_masked.where(SST_present_event_masked > 273)
        #     SST_future_event_masked = get_masked_ensemble(
        #         data_dict[f"{region}_future"], region, "future", "SST", event_mask_time, member_ids
        #     )
        #     SST_future_event_masked = SST_future_event_masked.where(SST_future_event_masked > 273)
        #     # Select the same extent as CAOI
        #     SST_present_event_masked = SST_present_event_masked.sel(
        #         lon=slice(extent[0], extent[1]), lat=slice(extent[2], extent[3])
        #     )
        #     SST_future_event_masked = SST_future_event_masked.sel(
        #         lon=slice(extent[0], extent[1]), lat=slice(extent[2], extent[3])
        #     )
        #     # Calculate the time mean for contour plotting
        #     SST_present_mean = calculate_time_mean(SST_present_event_masked)
        #     SST_future_mean = calculate_time_mean(SST_future_event_masked)
        #     # Plot contour lines (e.g., every 2K from 273 to 305)
        #     SST_contour_neg = ax.contour(
        #         SST_present_mean.lon, SST_present_mean.lat, 
        #         SST_future_mean - SST_present_mean,
        #         levels=[-4, -3, -2, -1],
        #         colors='black', linestyles='dashed', linewidths=2,
        #         transform=ccrs.PlateCarree()
        #     )
        #     SST_contour_pos = ax.contour(
        #         SST_present_mean.lon, SST_present_mean.lat, 
        #         SST_future_mean - SST_present_mean,
        #         levels=[1, 2, 3, 4, 5, 6],
        #         colors='black', linewidths=2,
        #         transform=ccrs.PlateCarree()
        #     )
        #     ax.clabel(SST_contour_neg, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
        #     ax.clabel(SST_contour_pos, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
        sea_ice_present_mean = calculate_time_mean(sea_ice_present)
        sea_ice_future_mean = calculate_time_mean(sea_ice_future)
        max_ice_concentration_present = sea_ice_present_mean.max().values
        max_ice_concentration_future = sea_ice_future_mean.max().values

        if max_ice_concentration_present >= 0.5:
            ax.contour(lons_plot, sea_ice_present_mean.lat, sea_ice_present_mean,
                    levels=[0.5], colors='green', linewidths=2,
                    transform=ccrs.PlateCarree(), alpha=1)
            print(f"Drew 50% sea ice contour for present (max: {max_ice_concentration_present:.3f})")
        else:
            print(f"No 50% sea ice contour available for present (max: {max_ice_concentration_present:.3f})")
        if max_ice_concentration_future >= 0.5:
            ax.contour(lons_plot, sea_ice_future_mean.lat, sea_ice_future_mean,
                    levels=[0.5], colors='green', linewidths=2,
                    transform=ccrs.PlateCarree(), alpha=1, linestyles='dashed')
            print(f"Drew 50% sea ice contour for future (max: {max_ice_concentration_future:.3f})")
        else:
            print(f"No 50% sea ice contour available for future (max: {max_ice_concentration_future:.3f})")

        # ax.contour(lons_SST, lats_SST, calculate_time_mean(SST_future), levels=[273], colors='green', linewidths=3, transform=ccrs.PlateCarree(), alpha=1, linestyles='dashed')
        # ax.contour(lons_SST, lats_SST, calculate_time_mean(SST_present), levels=[273], colors='green', linewidths=3, transform=ccrs.PlateCarree(), alpha=1)
        
        # if varname in ['PRECIP', 'EVAP', 'EP']:
        #     ax.contour(lons, lats, calculate_time_mean(SLP_future), levels=np.arange(-6, 8, 2), colors='black', linewidths=2, transform=ccrs.PlateCarree(), alpha=1, linestyles='dashed')
        #     ax.contour(lons, lats, calculate_time_mean(SLP_present), levels=np.arange(-6, 8, 2), colors='black', linewidths=2, transform=ccrs.PlateCarree(), alpha=1)
    if significance_mask is not None:
        # Create stippling where significant
        lon_grid, lat_grid = np.meshgrid(lons_plot, lats)
        
        # Only show stippling where significant (True values)
        significant_points = np.where(significance_mask.values)
        if len(significant_points[0]) > 0:  # Check if there are any significant points
            ax.scatter(lon_grid[significant_points], lat_grid[significant_points], 
                      s=3, c='black', marker='.', alpha=0.8,
                      transform=ccrs.PlateCarree(), label='p < 0.05')
            print(f"Added stippling for {len(significant_points[0])} significant points")
        else:
            print("No significant points found for stippling")   
    
    # cbar = plt.colorbar(im, ax=ax, orientation="horizontal", shrink=0.6, pad=0.05, extend=extend)
    divider = make_axes_locatable(ax)
    shrink = get_colorbar_shrink(extent)
    if extent == map_extent[region]:
        shrink = 0.9
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', shrink = shrink, pad = 0.05, extend=extend, alpha = 0.3)
    cbar.set_label(varname, size=15.5)
    if labname is not None:
        cbar.set_label(labname, size=15.5)
    cbar.ax.tick_params(labelsize=10.34)
    if cbarticks is not None:
        cbar.set_ticks(cbarticks)
    else:
        cbar.set_ticks(boundaries)
    # plt.title(title)

    if lag != 0:
        lag_name = f'_plus{lag}'
    else:
        lag_name = ''

    box_extent = box[f'{region}']
    # Limit data_2d, lons, and lats to the box extent
    if hasattr(data_2d, 'sel'):
        # If data_2d is an xarray DataArray
        data_2d_box = data_2d.sel(lon=slice(box_extent[0], box_extent[1]), lat=slice(box_extent[2], box_extent[3]))
    else:
        # If data_2d is a numpy array, mask lons/lats accordingly
        lon_mask = (lons >= box_extent[0]) & (lons <= box_extent[1])
        lat_mask = (lats >= box_extent[2]) & (lats <= box_extent[3])
        data_2d_box = data_2d[np.ix_(lat_mask, lon_mask)]
        data_2d_box = xr.DataArray(
            data_2d_box,
            coords={'lat': lats[lat_mask], 'lon': lons[lon_mask]},
            dims=['lat', 'lon']
        )
    print('Mean over box:', calculate_space_mean(data_2d_box))

    # if varname in ['CAOI']:
    #     # Plot rectangle of box extent
    #     ax.add_patch(plt.Rectangle(
    #         (box_extent[0], box_extent[2]),
    #         box_extent[1] - box_extent[0],
    #         box_extent[3] - box_extent[2],
    #         fill=False, edgecolor='black', linewidth=4,
    #         transform=ccrs.PlateCarree(), zorder=100
    #     ))

    if save_path == True:
        save_path = f"{fig_path}/{varname}_{timeframe_name}_{region}{lag_name}.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
        plt.savefig(save_path.replace('.png', '.pdf'), bbox_inches='tight')
        print(f"Saved plot to: {save_path}")



    plt.show()
    plt.close()

    return ax, cbar

def fdr(pvalues,alpha):
   """Given an array of p-values and a false discovery rate (FDR) level,
   return the p-value threshold that defines significance according to this FDR level.
   See Wilks (2016) for more info.
   """
   pvalues=pvalues[~np.isnan(pvalues)]
   pvalues=pvalues[~np.isinf(pvalues)]

   sortidx = np.argsort(pvalues)
   psorted = pvalues[sortidx]
   nval = len(pvalues)
   ifdr = np.argmax((psorted < alpha*np.arange(1,nval+1)/nval)[::-1])
   if ifdr == 0 and psorted[-1]>= alpha:
       ifdr=nval-1
   ifdr = nval - ifdr - 1
   return pvalues[sortidx[ifdr]]

def fdr_xr_2d(pvalues,alpha):
   """Given an map of p-values and a false discovery rate (FDR) level,
   return the mask that defines significance according to this FDR level.
   See Wilks (2016) for more info.
   """    
   lats = pvalues.lat
   lons = pvalues.lon
   pvalues=np.array(pvalues)
   pthresh = fdr(pvalues.reshape(-1),alpha)
   significance_mask_da = xr.DataArray(
        pvalues<pthresh,
        coords={'lat': lats, 'lon': lons},
        dims=['lat', 'lon'], 
        name='significance'
    )
   return significance_mask_da

def compute_significance_ttest_map(present_data, future_data, alpha=0.05):
    """
    Compute statistical significance between present and future data using t-test.
    
    Parameters:
    -----------
    present_data : xarray.DataArray
        Present data with dimensions (time, lat, lon)
    future_data : xarray.DataArray 
        Future data with dimensions (time, lat, lon)
    alpha : float
        Significance level (default 0.05)
        
    Returns:
    --------
    significance_mask : xarray.DataArray
        Boolean mask where True indicates significant differences
    p_values : xarray.DataArray
        P-values from t-tests at each grid point
    """
    # Get dimensions
    lats = present_data.lat
    lons = present_data.lon
    
    # Initialize output arrays
    p_values = np.full((len(lats), len(lons)), np.nan)
    significance_mask = np.full((len(lats), len(lons)), False, dtype=bool)
    
    # Perform t-test at each grid point
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            present_vals = present_data.sel(lat=lat, lon=lon).values
            future_vals = future_data.sel(lat=lat, lon=lon).values
            
            # Remove NaN values
            present_clean = present_vals[~np.isnan(present_vals)]
            future_clean = future_vals[~np.isnan(future_vals)]
            
            # Only perform test if we have enough data
            if len(present_clean) > 3 and len(future_clean) > 3:
                try:
                    # Two-sample t-test (assuming unequal variances)
                    t_stat, p_val = ttest_ind(future_clean, present_clean, equal_var=False)
                    p_values[i, j] = p_val
                    significance_mask[i, j] = p_val < alpha
                except:
                    # Handle any statistical errors
                    p_values[i, j] = np.nan
                    significance_mask[i, j] = False
    
    # Convert to xarray DataArrays
    p_values_da = xr.DataArray(
        p_values, 
        coords={'lat': lats, 'lon': lons}, 
        dims=['lat', 'lon'],
        name='p_values'
    )
    
    significance_mask_da = xr.DataArray(
        significance_mask,
        coords={'lat': lats, 'lon': lons},
        dims=['lat', 'lon'], 
        name='significance'
    )
    
    return significance_mask_da, p_values_da

def plot_map_10th_with_significance(region, member_ids=member_ids, varname='CAOI', lag=0, alpha=0.05):
    """
    Enhanced version of plot_map_10th with built-in significance testing and stippling.
    """
    # Get the configuration for this variable
    if varname not in VAR_CONFIG:
        raise ValueError(f"Variable '{varname}' not found in VAR_CONFIG. Available variables: {list(VAR_CONFIG.keys())}")
    
    if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
        data_dict_v = data_dict_precip
        if lag == 0:
            event_mask_time_func = event_plus0_mask_time
        elif lag == 1:
            event_mask_time_func = event_plus1_mask_time
        elif lag == 2:
            event_mask_time_func = event_plus2_mask_time
        elif lag == 3:
            event_mask_time_func = djf_mask
    else:
        data_dict_v = data_dict
        event_mask_time_func = event_mask_time

    config = VAR_CONFIG[varname]
    diff_config = VAR_CONFIG_DIFF[varname]

    # Get masked ensemble data
    present_data = get_masked_ensemble(data_dict_v[f'{region}_present'], region, "present", varname, event_mask_time_func, member_ids)
    future_data = get_masked_ensemble(data_dict_v[f'{region}_future'], region, "future", varname, event_mask_time_func, member_ids)

    # Debug the data first
    # print("Debugging data before significance testing...")
    # debug_significance_data(present_data, future_data)

    # Calculate time means
    present_mean_time = calculate_time_mean(present_data)
    future_mean_time = calculate_time_mean(future_data)
    diff_mean_time = future_mean_time - present_mean_time

    # Perform significance testing
    print(f"Performing significance testing for {varname} in {region}...")
    significance_mask, p_values = compute_significance_ttest_map(present_data, future_data, alpha)
    significance_mask = fdr_xr_2d(p_values, alpha=0.1)  # Apply FDR correction
    
    print(f"Significance testing results:")
    print(f"Total grid points: {significance_mask.size}")
    print(f"Significant points (p < {alpha}): {significance_mask.sum().item()}")
    print(f"Percentage significant: {(significance_mask.sum() / significance_mask.size * 100).item():.1f}%")

    # Set up plotting coordinates
    lonmin, lonmax = present_data.lon.min(), present_data.lon.max()
    latmin, latmax = present_data.lat.min(), present_data.lat.max()
    extent = [lonmin, lonmax, latmin, latmax]
    lons = present_data.lon
    lats = present_data.lat

    # Plot present (no stippling)
    ax_present = plot_single_map(
        data_2d=present_mean_time,
        lons=lons, lats=lats, region=region,
        timeframe_name='10th_present', extent=extent,
        cmap_name=config['cmap_name'], boundaries=config['boundaries'],
        cbarticks=config['cbarticks'], varname=varname,
        labname=config['labname'], extend='max',
        save_path=True, lag=lag,
        significance_mask=None  # No stippling for present
    )

    # Plot future (no stippling)
    ax_future = plot_single_map(
        data_2d=future_mean_time,
        lons=lons, lats=lats, region=region,
        timeframe_name='10th_future', extent=extent,
        cmap_name=config['cmap_name'], boundaries=config['boundaries'],
        cbarticks=config['cbarticks'], varname=varname,
        labname=config['labname'], extend='max',
        save_path=True, lag=lag,
        significance_mask=None  # No stippling for future
    )

    # Plot difference WITH stippling
    ax_diff = plot_single_map(
        data_2d=diff_mean_time,
        lons=lons, lats=lats, region=region,
        timeframe_name='10th_diff', extent=extent,
        cmap_name=diff_config['cmap_name'], boundaries=diff_config['boundaries'],
        cbarticks=diff_config['cbarticks'], varname=varname,
        labname=diff_config['labname'], extend='both',
        save_path=True, lag=lag,
        significance_mask=significance_mask  # Add stippling for differences
    )

    return significance_mask, p_values, ax_present, ax_future, ax_diff

def compute_significance_frequency(present_data, future_data, threshold=8, alpha=0.05):
    """
    Compute statistical significance for frequency differences using z-test for proportions.
    
    Parameters:
    -----------
    present_data : xarray.DataArray
        Present data with dimensions (time, lat, lon)
    future_data : xarray.DataArray 
        Future data with dimensions (time, lat, lon)
    threshold : float
        Threshold value (default 8 for CAOI > 8K)
    alpha : float
        Significance level (default 0.05)
        
    Returns:
    --------
    significance_mask : xarray.DataArray
        Boolean mask where True indicates significant differences
    p_values : xarray.DataArray
        P-values from z-tests at each grid point
    """
    
    # Get dimensions
    lats = present_data.lat
    lons = present_data.lon
    
    # Initialize output arrays
    p_values = np.full((len(lats), len(lons)), np.nan)
    significance_mask = np.full((len(lats), len(lons)), False, dtype=bool)
    
    # Get total counts
    n_present = len(present_data.time)
    n_future = len(future_data.time)
    
    print(f"Sample sizes: Present={n_present}, Future={n_future}")
    
    # Perform z-test at each grid point
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            present_vals = present_data.sel(lat=lat, lon=lon).values
            future_vals = future_data.sel(lat=lat, lon=lon).values
            
            # Remove NaN values and count exceedances
            present_clean = present_vals[~np.isnan(present_vals)]
            future_clean = future_vals[~np.isnan(future_vals)]
            
            if len(present_clean) > 5 and len(future_clean) > 5:
                # Count exceedances
                x1 = np.sum(present_clean > threshold)  # present exceedances
                x2 = np.sum(future_clean > threshold)   # future exceedances
                n1 = len(present_clean)  # present sample size
                n2 = len(future_clean)   # future sample size
                
                # Calculate proportions
                p1 = x1 / n1
                p2 = x2 / n2
                
                # Pooled proportion
                p_pool = (x1 + x2) / (n1 + n2)
                
                # Only test if there are some exceedances
                if x1 + x2 > 0:
                    try:
                        # Standard error for difference in proportions
                        se = np.sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
                        
                        if se > 0:
                            # Z-statistic
                            z = (p2 - p1) / se
                            
                            # Two-tailed p-value
                            p_val = 2 * (1 - norm.cdf(np.abs(z)))
                            
                            p_values[i, j] = p_val
                            significance_mask[i, j] = p_val < alpha
                        else:
                            p_values[i, j] = np.nan
                            significance_mask[i, j] = False
                    except:
                        p_values[i, j] = np.nan
                        significance_mask[i, j] = False
                else:
                    # No exceedances in either period
                    p_values[i, j] = 1.0
                    significance_mask[i, j] = False
    
    # Convert to xarray DataArrays
    p_values_da = xr.DataArray(
        p_values, 
        coords={'lat': lats, 'lon': lons}, 
        dims=['lat', 'lon'],
        name='p_values_frequency'
    )
    
    significance_mask_da = xr.DataArray(
        significance_mask,
        coords={'lat': lats, 'lon': lons},
        dims=['lat', 'lon'], 
        name='significance_frequency'
    )
    
    return significance_mask_da, p_values_da

def plot_map_frequency_10th_with_significance(region, alpha=0.05):
    """
    Enhanced version of plot_map_frequency_10th with significance testing and stippling.
    """
    CAOI_present = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', event_mask_time, member_ids)
    CAOI_future = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', event_mask_time, member_ids)

    # Debug the data
    print(f"CAOI present shape: {CAOI_present.shape}")
    print(f"CAOI future shape: {CAOI_future.shape}")
    print(f"CAOI present time steps: {len(CAOI_present.time)}")
    print(f"CAOI future time steps: {len(CAOI_future.time)}")

    # Perform significance testing
    print(f"Performing significance testing for frequency maps in {region}...")
    significance_mask, p_values = compute_significance_frequency(CAOI_present, CAOI_future, threshold=8, alpha=alpha)
    significance_mask = fdr_xr_2d(p_values, alpha=0.1)  # Apply FDR correction
    
    print(f"Significance testing results:")
    print(f"Total grid points: {significance_mask.size}")
    print(f"Significant points (p < {alpha}): {significance_mask.sum().item()}")
    print(f"Percentage significant: {(significance_mask.sum() / significance_mask.size * 100).item():.1f}%")

    lonmin, lonmax = CAOI_present.lon.min(), CAOI_present.lon.max()
    latmin, latmax = CAOI_present.lat.min(), CAOI_present.lat.max()
    box = [lonmin, lonmax, latmin, latmax]
    lons = CAOI_present.lon
    lats = CAOI_present.lat

    CAOI_present_nans = xr.DataArray(
        np.all(np.isnan(CAOI_present), axis=0),
        coords={'lat': lats, 'lon': lons},
        dims=['lat', 'lon']
    )
    CAOI_future_nans = xr.DataArray(
        np.all(np.isnan(CAOI_future), axis=0),
        coords={'lat': lats, 'lon': lons},
        dims=['lat', 'lon']
    )

    # Calculate frequencies
    CAOI_present_count = np.nansum(CAOI_present > 8, axis=(0))
    CAOI_future_count = np.nansum(CAOI_future > 8, axis=(0))

    CAOI_present_freq = CAOI_present_count / len(CAOI_present.time) * 100
    CAOI_future_freq = CAOI_future_count / len(CAOI_future.time) * 100
    CAOI_diff_freq = CAOI_future_freq - CAOI_present_freq

    print('max present', np.max(CAOI_present_freq))
    print('mean present', np.nanmean(CAOI_present_freq))
    print('max future', np.max(CAOI_future_freq))
    print('mean future', np.nanmean(CAOI_future_freq))
    print('max diff', np.max(CAOI_diff_freq))
    print('mean diff', np.nanmean(CAOI_diff_freq))

    CAOI_present_freq = np.where(CAOI_present_freq == 0, np.nan, CAOI_present_freq)
    CAOI_future_freq = np.where(CAOI_future_freq == 0, np.nan, CAOI_future_freq)
    CAOI_diff_freq = np.where(CAOI_diff_freq == 0, np.nan, CAOI_diff_freq)
    CAOI_diff_freq = np.where(CAOI_present_nans.values | CAOI_future_nans.values, np.nan, CAOI_diff_freq)

    # Plot present (stipple grid cells where CAOI_present is always NaN)
    plot_single_map(
        data_2d=CAOI_present_freq,
        lons=lons, lats=lats, region=region,
        timeframe_name='present', extent=box,
        cmap_name='bilbao', boundaries=np.arange(0, 110, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='frequency exceedance 8K [%]',
        extend='neither', save_path=True,
        significance_mask=CAOI_present_nans
    )

    # Plot future (stipple grid cells where CAOI_future is always NaN)
    plot_single_map(
        data_2d=CAOI_future_freq,
        lons=lons, lats=lats, region=region,
        timeframe_name='future', extent=box,
        cmap_name='bilbao', boundaries=np.arange(0, 110, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='frequency exceedance 8K [%]',
        extend='neither', save_path=True,
        significance_mask=CAOI_future_nans
    )

    # Plot difference WITH stippling
    plot_single_map(
        data_2d=CAOI_diff_freq,
        lons=lons, lats=lats, region=region,
        timeframe_name='diff', extent=box,
        cmap_name='vik', boundaries=np.arange(-30, 35, 5),
        cbarticks=np.arange(-30, 40, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='Δ frequency exceedance 8K [%]',
        extend='both', save_path=True,
        significance_mask=significance_mask  # Add stippling for differences
    )

    return significance_mask, p_values

def compute_significance_ttest(present_data, future_data, alpha=0.05):
    """
    Compute t-test significance between present and future data.
    
    Parameters:
    -----------
    present_data : array-like
        Present data values
    future_data : array-like
        Future data values
    alpha : float
        Significance level (default 0.05)
        
    Returns:
    --------
    is_significant : bool
        True if difference is significant
    p_value : float
        P-value from t-test
    t_statistic : float
        T-statistic from t-test
    """
    # Remove NaN values
    present_clean = present_data[~np.isnan(present_data)]
    future_clean = future_data[~np.isnan(future_data)]
    
    # Only perform test if we have enough data
    if len(present_clean) > 3 and len(future_clean) > 3:
        try:
            # Two-sample t-test (assuming unequal variances)
            t_stat, p_val = ttest_ind(future_clean, present_clean, equal_var=False)
            return p_val < alpha, p_val, t_stat
        except:
            return False, np.nan, np.nan
    else:
        return False, np.nan, np.nan

def plot_caoi_contribution_combined_with_significance(region):
    # Enforce a concrete sans-serif font for this figure only.
    # In interactive sessions, Matplotlib may keep a stale font cache, so add Arial explicitly if available.
    arial_user_font = os.path.expanduser("~/.local/share/fonts/truetype/arial/arial.ttf")
    if os.path.exists(arial_user_font):
        try:
            matplotlib.font_manager.fontManager.addfont(arial_user_font)
        except Exception:
            pass
    available_fonts = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
    preferred_font = "Arial" if "Arial" in available_fonts else "Liberation Sans"

    if region == "JS":
        ylims_temp = (245, 305)
        ylims_caoi = (-5, 25)
    elif region == "NA":
        ylims_temp = (255, 310)
        ylims_caoi = (-15, 40)

    ylims_temp = (270, 310)
    ylims_caoi = (-11, 16)

    # EVENT DATA (foreground)
    event_CAOI_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    event_CAOI_present = calculate_space_mean(event_CAOI_present)
    event_SSTH_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    event_SSTH_present_space_mean = calculate_space_mean(event_SSTH_present)
    event_SSTH_present_space_mean_noweights = calculate_space_mean_no_weights(event_SSTH_present)

    print(f"Event SSTH present space mean (weighted): {event_SSTH_present_space_mean.values.flatten().mean():.2f} K")
    print(f"Event SSTH present space mean (unweighted): {event_SSTH_present_space_mean_noweights.values.flatten().mean():.2f} K")


    event_SSTH_present = calculate_space_mean(event_SSTH_present)
    event_TH_850_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)
    event_TH_850_present = calculate_space_mean(event_TH_850_present)

    event_CAOI_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    event_CAOI_future = calculate_space_mean(event_CAOI_future)
    event_SSTH_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    event_SSTH_future = calculate_space_mean(event_SSTH_future)
    event_TH_850_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)
    event_TH_850_future = calculate_space_mean(event_TH_850_future)

    # DJF DATA (background)
    djf_CAOI_present = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    djf_CAOI_present = calculate_space_mean(djf_CAOI_present)
    djf_SSTH_present = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    djf_SSTH_present = calculate_space_mean(djf_SSTH_present)
    djf_TH_850_present = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)
    djf_TH_850_present = calculate_space_mean(djf_TH_850_present)

    djf_CAOI_future = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    djf_CAOI_future = calculate_space_mean(djf_CAOI_future)
    djf_SSTH_future = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    djf_SSTH_future = calculate_space_mean(djf_SSTH_future)
    djf_TH_850_future = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)
    djf_TH_850_future = calculate_space_mean(djf_TH_850_future)

    # Prepare data for boxplots - EVENT DATA
    event_data_temp = [
        event_SSTH_present.values.flatten(),
        event_SSTH_future.values.flatten(),
        event_TH_850_present.values.flatten(),
        event_TH_850_future.values.flatten(),
    ]
    event_data_caoi = [
        event_CAOI_present.values.flatten(),
        event_CAOI_future.values.flatten(),
    ]

    # Prepare data for boxplots - DJF DATA
    djf_data_temp = [
        djf_SSTH_present.values.flatten(),
        djf_SSTH_future.values.flatten(),
        djf_TH_850_present.values.flatten(),
        djf_TH_850_future.values.flatten(),
    ]
    djf_data_caoi = [
        djf_CAOI_present.values.flatten(),
        djf_CAOI_future.values.flatten(),
    ]

    # Filter out NaN values
    event_data_temp = [d[~np.isnan(d)] for d in event_data_temp]
    event_data_caoi = [d[~np.isnan(d)] for d in event_data_caoi]
    djf_data_temp = [d[~np.isnan(d)] for d in djf_data_temp]
    djf_data_caoi = [d[~np.isnan(d)] for d in djf_data_caoi]

    # Calculate mean differences (future - present) for EVENT data
    mean_diff_SSTH = np.nanmean(event_data_temp[1]) - np.nanmean(event_data_temp[0])
    mean_diff_TH850 = np.nanmean(event_data_temp[3]) - np.nanmean(event_data_temp[2])
    mean_diff_CAOI = np.nanmean(event_data_caoi[1]) - np.nanmean(event_data_caoi[0])

    # Mean diffs for DJF data
    mean_diff_SSTH_DJF = np.nanmean(djf_data_temp[1]) - np.nanmean(djf_data_temp[0])
    mean_diff_TH850_DJF = np.nanmean(djf_data_temp[3]) - np.nanmean(djf_data_temp[2])
    mean_diff_CAOI_DJF = np.nanmean(djf_data_caoi[1]) - np.nanmean(djf_data_caoi[0])

    # Perform significance testing for EVENT data
    is_sig_SSTH, p_val_SSTH, t_stat_SSTH = compute_significance_ttest(event_data_temp[0], event_data_temp[1])
    is_sig_TH850, p_val_TH850, t_stat_TH850 = compute_significance_ttest(event_data_temp[2], event_data_temp[3])
    is_sig_CAOI, p_val_CAOI, t_stat_CAOI = compute_significance_ttest(event_data_caoi[0], event_data_caoi[1])

    # Print significance results
    print(f"Significance testing results for {region}:")
    print(f"SSTH: Δ={mean_diff_SSTH:.2f}K, p={p_val_SSTH:.4f}, significant={is_sig_SSTH}")
    print(f"TH850: Δ={mean_diff_TH850:.2f}K, p={p_val_TH850:.4f}, significant={is_sig_TH850}")
    print(f"CAOI: Δ={mean_diff_CAOI:.2f}K, p={p_val_CAOI:.4f}, significant={is_sig_CAOI}")
    print("--------------------------------------------------")
    print(compute_significance_ttest(event_data_temp[0], event_data_temp[1]))
    print(compute_significance_ttest(event_data_temp[2], event_data_temp[3]))
    print(compute_significance_ttest(event_data_caoi[0], event_data_caoi[1]))
    print("==================================================")
    print(f"DJF SSTH: Δ={mean_diff_SSTH_DJF:.2f}K")
    print(f"DJF TH850: Δ={mean_diff_TH850_DJF:.2f}K")
    print(f"DJF CAOI: Δ={mean_diff_CAOI_DJF:.2f}K")
    print("==================================================")


    # Create labels with asterisks for significant differences
    asterisk_SSTH = "*" if is_sig_SSTH else ""
    asterisk_TH850 = "*" if is_sig_TH850 else ""
    asterisk_CAOI = "*" if is_sig_CAOI else ""

    # --- Plot both boxplots in the same figure ---
    with plt.rc_context({
        "font.family": "sans-serif",
        "font.sans-serif": [preferred_font, "Arial", "Liberation Sans", "DejaVu Sans"],
    }):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 9), gridspec_kw={'width_ratios': [2, 1]})

        # SSTH & TH850 - Plot DJF first (background with alpha=0.2)
        bp_temp_djf = ax1.boxplot(
            djf_data_temp,
            positions=[0, 1, 2, 3],
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            showmeans=True,
            meanline=True,
            meanprops={'color': 'black', 'linewidth': 1.5, 'linestyle': '--', 'alpha': 0.2},
            medianprops={'linewidth': 0, 'alpha': 0.2},
        )
        colors_temp_djf = ["lightblue", "orange", "lightblue", "orange"]
        for patch, color in zip(bp_temp_djf["boxes"], colors_temp_djf):
            patch.set_facecolor(color)
            patch.set_alpha(0.2)
        
        # Set line properties for DJF boxplot to be less prominent
        for element in ['whiskers', 'fliers', 'medians', 'caps']:
            for line in bp_temp_djf[element]:
                line.set_alpha(0.2)
        
        # SSTH & TH850 - Plot events on top (foreground)
        bp_temp_event = ax1.boxplot(
            event_data_temp,
            positions=[0, 1, 2, 3],
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            showmeans=True,
            meanline=True,
            meanprops={'color': 'black', 'linewidth': 1.5, 'linestyle': '--'},
            medianprops={'linewidth': 0},
        )
        colors_temp_event = ["blue", "orange", "blue", "orange"]
        for patch, color in zip(bp_temp_event["boxes"], colors_temp_event):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax1.set_ylabel("Θ [K]", fontsize=20, fontname=preferred_font)
        ax1.set_ylim(ylims_temp)
        ax1.tick_params(axis="both", labelsize=16)
        ax1.grid(axis="y", linestyle="--", alpha=0.5)
        ax1.set_xticks([0.5, 2.5])

        SST_text = 'SST'
        TH_text = '850'
        
        # Update labels to include asterisks for significant differences
        xticklabels = [
            f"$\\Theta_{{{SST_text}}}$\nΔ={mean_diff_SSTH:.2f}K{asterisk_SSTH}",
            f"$\\Theta_{{{TH_text}}}$\nΔ={mean_diff_TH850:.2f}K{asterisk_TH850}"
        ]
        ax1.set_xticklabels(xticklabels, fontsize=20, fontname=preferred_font)
        print(ax1.get_ylim())

        # CAOI - Plot DJF first (background with alpha=0.2)
        bp_caoi_djf = ax2.boxplot(
            djf_data_caoi,
            positions=[0, 1],
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            showmeans=True,
            meanline=True,
            meanprops={'color': 'black', 'linewidth': 1.5, 'linestyle': '--', 'alpha': 0.2},
            medianprops={'linewidth': 0, 'alpha': 0.2},
        )
        colors_caoi_djf = ["lightblue", "orange"]
        for patch, color in zip(bp_caoi_djf["boxes"], colors_caoi_djf):
            patch.set_facecolor(color)
            patch.set_alpha(0.2)
        
        # Set line properties for DJF boxplot to be less prominent
        for element in ['whiskers', 'fliers', 'medians', 'caps']:
            for line in bp_caoi_djf[element]:
                line.set_alpha(0.2)

        # CAOI - Plot events on top (foreground)
        bp_caoi_event = ax2.boxplot(
            event_data_caoi,
            positions=[0, 1],
            patch_artist=True,
            showfliers=False,
            widths=0.6,
            showmeans=True,
            meanline=True,
            meanprops={'color': 'black', 'linewidth': 1.5, 'linestyle': '--'},
            medianprops={'linewidth': 0},
        )
        colors_caoi_event = ["blue", "orange"]
        for patch, color in zip(bp_caoi_event["boxes"], colors_caoi_event):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Move ylabel to the right side
        ax2.set_ylabel("CAOI [K]", fontsize=20, labelpad=10, fontname=preferred_font)
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.set_ylim(ylims_caoi)
        ax2.tick_params(axis="y", labelsize=16)
        ax2.grid(axis="y", linestyle="--", alpha=0.5)
        ax2.set_xticks([0.5])
        
        # Update CAOI label to include asterisk for significant difference
        ax2.set_xticklabels([f"CAOI\nΔ={mean_diff_CAOI:.2f}K{asterisk_CAOI}"], fontsize=20, fontname=preferred_font)
        print(ax2.get_ylim())

        # Force tick/legend fonts for all text artists in this figure.
        for ax in (ax1, ax2):
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontname(preferred_font)

        # Add legend to distinguish between DJF and events
        legend_elements = [
            # Patch(facecolor='gray', alpha=0.2, label='DJF'),
            # Patch(facecolor='gray', alpha=0.7, label='events'),
            Patch(facecolor='blue', alpha=0.7, label='present'),
            Patch(facecolor='orange', alpha=0.7, label='future')
        ]
        legend = ax1.legend(handles=legend_elements, loc='upper left')
        for text in legend.get_texts():
            text.set_fontname(preferred_font)

        # # Debug: show the actual font files Matplotlib resolves for this figure's text.
        # y1_props = ax1.yaxis.label.get_fontproperties()
        # y2_props = ax2.yaxis.label.get_fontproperties()
        # legend_props = legend.get_texts()[0].get_fontproperties()
        # print(f"[FONT DEBUG] preferred_font={preferred_font}")
        # print(f"[FONT DEBUG] ax1 ylabel family={y1_props.get_family()} name={y1_props.get_name()}")
        # print(f"[FONT DEBUG] ax1 ylabel file={matplotlib.font_manager.findfont(y1_props, fallback_to_default=True)}")
        # print(f"[FONT DEBUG] ax2 ylabel family={y2_props.get_family()} name={y2_props.get_name()}")
        # print(f"[FONT DEBUG] ax2 ylabel file={matplotlib.font_manager.findfont(y2_props, fallback_to_default=True)}")
        # print(f"[FONT DEBUG] legend text family={legend_props.get_family()} name={legend_props.get_name()}")
        # print(f"[FONT DEBUG] legend text file={matplotlib.font_manager.findfont(legend_props, fallback_to_default=True)}")
        # ticklabels = ax1.get_xticklabels()
        # if ticklabels:
        #     tick_props = ticklabels[0].get_fontproperties()
        #     print(f"[FONT DEBUG] ax1 xtick[0] family={tick_props.get_family()} name={tick_props.get_name()}")
        #     print(f"[FONT DEBUG] ax1 xtick[0] file={matplotlib.font_manager.findfont(tick_props, fallback_to_default=True)}")
        # yticklabels = ax1.get_yticklabels()
        # if yticklabels:
        #     ytick_props = yticklabels[0].get_fontproperties()
        #     print(f"[FONT DEBUG] ax1 ytick[0] family={ytick_props.get_family()} name={ytick_props.get_name()}")
        #     print(f"[FONT DEBUG] ax1 ytick[0] file={matplotlib.font_manager.findfont(ytick_props, fallback_to_default=True)}")

        # # Add note about significance testing
        # fig.suptitle(f'{region} - * indicates p < 0.05', fontsize=14, y=0.02)

        plt.tight_layout()
        plt.savefig(f"{fig_path}/CAOI_contribution_combined_overlay_{region}_with_significance.png", 
                    dpi=300, bbox_inches="tight", facecolor="white", transparent=False)
        plt.savefig(f"{fig_path}/CAOI_contribution_combined_overlay_{region}_with_significance.pdf", bbox_inches='tight')
        plt.show()


class LazyEnsemble(dict):
    """Ensemble data of one region and timeframe (key "JS_present", "NA_future", ...),
    read from disk the first time it is used."""

    def __init__(self, reader):
        super().__init__()
        self.reader = reader

    def __missing__(self, key):
        region, timeframe = key.split("_")
        self[key] = self.reader(region, timeframe, member_ids)
        return self[key]


# Ensemble data read by the plotting functions above
data_dict = LazyEnsemble(read_ensemble_data)
data_dict_precip = LazyEnsemble(read_ensemble_data_precip)
data_dict_ice = LazyEnsemble(read_ensemble_data_ice)


def free_memory():
    """Release the ensemble data after a figure to keep the memory use low."""
    for data in (data_dict, data_dict_precip, data_dict_ice):
        data.clear()
    gc.collect()


def maps(regions):
    """Maps of CAOI intensity, CAO frequency and evaporation (present, future, difference)."""
    for region in regions:
        plot_map_10th_with_significance(region, varname='CAOI')
        free_memory()
        plot_map_frequency_10th_with_significance(region)
        free_memory()
        plot_map_10th_with_significance(region, varname='EVAP')
        free_memory()


def contribution(regions):
    """Contribution of the CAOI components to the future change."""
    for region in regions:
        plot_caoi_contribution_combined_with_significance(region)
        free_memory()


def contribution_JS_northern_box_50N():
    """As the JS contribution figure, but with the northern box boundary at 50N."""
    global fig_path
    original_box, original_path = box["JS"], fig_path
    box["JS"] = [original_box[0], original_box[1], original_box[2], 50]
    fig_path = f"{original_path}/S3"
    os.makedirs(fig_path, exist_ok=True)
    try:
        plot_caoi_contribution_combined_with_significance("JS")
        free_memory()
    finally:
        box["JS"], fig_path = original_box, original_path


# Figure name -> function that produces it
FIGURES = {
    "F1": lambda: maps(["JS", "NA"]),         # panels with the future-present difference
    "S1": lambda: maps(["JS"]),               # present and future panels, Japan Sea
    "S2": lambda: maps(["NA"]),               # present and future panels, northwest Atlantic
    "F2": lambda: contribution(["JS", "NA"]),
    "S3": contribution_JS_northern_box_50N,
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

# %%
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from lagranto.plotting import plot_trajs
from lagranto import Tra 
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
import matplotlib.ticker as mticker
import matplotlib.path as mpath
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm
import matplotlib
import pandas as pd
import netCDF4 as nc
import numpy as np
import datetime as dt
from scipy.stats import norm
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from scipy.stats import ttest_ind

from matplotlib.lines import Line2D  # Make sure this import is included
import xarray as xr
import matplotlib.dates as mdates  # Import mdates for date formatting
import cftime
import os
from cartopy.io import shapereader
from cartopy.feature import LAND
import cftime  # Add this import
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable


# %%
member_ids = [1000, 1100, 1200, 1300, 1400]
regions = ["JS", "NA"]
timeframes = ["present", "future"]
base_path = f"/net/litho/atmosdyn2/nfieldhouse/data/cao/clim/"
fig_path = '/home/nfieldhouse/ma/plots/statistics/ensemble'

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
map_extent = {
    "JS": [120, 160, 20.2617801, 69.26701571],  # for map
    "NA": [275, 305, 20.2617801, 54.18848168]   # for map
}

map_extent = {
    "JS": [120, 160, 20.2617801, 64.554974],  # for map
    "NA": [270, 310, 20.2617801, 64.554974]   # for map
}

map_extent_precip = {
    "JS": [120+10, 260-10, 20+5, 70],  # for map, central lon = 190
    "NA": [275, 5, 20, 75]   # for map, central lon = -40 (320)
}

map_extent_precip = {
    'JS': [130, 190, 25, 70], # central lon = 160
    'NA': [275, 335, 25, 70]  # central lon = -55 (305)
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

# Add this configuration dictionary at the top of your file or before the function
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
        'labname': '$\Theta_{SST}$ [K]'
    },
    'TH_850_masked': {
        'cmap_name': 'plasma',
        'boundaries': np.arange(250, 290, 5),
        'cbarticks': np.arange(250, 300, 10),
        'labname': '$\Theta_{850}$ [K]'
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
        'labname': '$\Theta_{SST}$ anomaly [K]'
    },
    'TH_850_ano': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(-15, 16, 1.5),
        'cbarticks': np.arange(-15, 16, 3),
        'labname': '$\Theta_{850}$ anomaly [K]'
    }
}

# Add difference configurations for the diff plots
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
        'labname': 'Δ $\Theta_{SST}$ [K]'
    },
    'TH_850_masked': {
        'cmap_name': 'RdBu_r',
        'boundaries': np.arange(-2, 2.5, 0.5),
        'cbarticks': np.arange(-2, 3, 1),
        'labname': 'Δ $\Theta_{850}$ [K]'
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
        'labname': 'Δ $\Theta_{SST}$ anomaly [K]'
    },
    'SSTH_ano': {
        'cmap_name': 'coolwarm',
        'boundaries': np.arange(-2, 2.2, 0.2),
        'cbarticks': np.arange(-2, 2.2, 0.4),
        'labname': '$\Theta_{850}$ anomaly [K]'
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

# Read data
def read_masked_CAOI_data(region, timeframe, member_ids):
    """
    Read CAOI data masked by sea-ice, land, box (time, lat, lon) for all ensemble members.
    Returns: CAOI_djf_ens, CAOI_djf_ano_ens with dims (time, lat, lon, member)
    """
    CAOI_djf_list = []
    CAOI_djf_ano_list = []
    for member_id in member_ids:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/"
        ds_C_ano = xr.open_dataset(path + f"yearly_hour_anomalies_{timeframe}", decode_times=True)
        ds_C = xr.open_dataset(path + f"all_times_{timeframe}", decode_times=True)

        CAOI = -ds_C.CAOI[:, 0, :, :]
        CAOI_djf = -ds_C.sel(time=djf_mask(ds_C, timeframe)).CAOI[:, 0, :, :]
        CAOI_ano = ds_C_ano.CAOI[:, 0, :, :]
        CAOI_djf_ano = ds_C_ano.sel(time=djf_mask(ds_C_ano, timeframe)).CAOI[:, 0, :, :]

        # Add member dimension
        CAOI_djf = CAOI_djf.expand_dims(member=[member_id])
        CAOI_djf_ano = CAOI_djf_ano.expand_dims(member=[member_id])

        CAOI_djf_list.append(CAOI_djf)
        CAOI_djf_ano_list.append(CAOI_djf_ano)

    # Concatenate along member dimension
    CAOI_djf_ens = xr.concat(CAOI_djf_list, dim='member')
    CAOI_djf_ano_ens = xr.concat(CAOI_djf_ano_list, dim='member')
    return CAOI_djf_ens, CAOI_djf_ano_ens

def write_data(region, timeframe, member_id):
    """Write DJF netcdf data masked by sea-ice, land, NOT BOX (time, lat, lon)"""
    if timeframe == "present":
        years = range(1990, 2001)  # for example: 1990 to 2000
        cesm_type = 'HIST'
    elif timeframe == "future":
        years = range(2090, 2101)
        cesm_type = 'SSP370'

    cesm = f"/net/meso/climphys/cesm212/b.e212.B{cesm_type}cmip6.f09_g17.{member_id}/archive/atm/hist/"

    print('Opening files')
    # Construct file paths for each year
    files_to_open = [
        os.path.join(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.{member_id}.cam.h3.{year}-01-01-00000.nc")
        for year in years
    ]
    print('Opened files')

    # Open and combine the datasets
    ds = xr.open_mfdataset(files_to_open, combine='by_coords', decode_times=True)
    SST = ds.sel(time=djf_mask(ds, timeframe), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).SST
    print('Opened SST')
    T = ds.sel(time=djf_mask(ds, timeframe), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).T[:,25,:,:]
    print('Opened T')
    PS = ds.sel(time=djf_mask(ds, timeframe), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).PS
    print('Opened PS')
    SHFLX = ds.sel(time=djf_mask(ds, timeframe), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).SHFLX
    print('Opened SHFLX')
    LHFLX = ds.sel(time=djf_mask(ds, timeframe), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).LHFLX
    print('Opened LHFLX')
    SSTH = SST*(1000/(PS/100))**(0.286)
    TH_850 = T*(1000/(85953/100))**(0.286)

    # exclude grid points with sea ice  (SST < 273K) & land
    sea_ice = np.where(SST < 273, np.nan,1)
    LSM = np.where(SST < 271,np.nan,1)
    mask = LSM * sea_ice
    CAOI = SSTH - TH_850

    CAOI_masked = mask*CAOI
    SSTH_masked = mask*SSTH
    TH_850_masked = mask*TH_850
    SHFLX_masked = mask*SHFLX
    LHFLX_masked = mask*LHFLX

    ds_to_save = xr.Dataset(
        {
            "CAOI": (["time", "lat", "lon"], CAOI_masked.data),
            "SST": (["time", "lat", "lon"], SST.data),
            "SSTH_masked": (["time", "lat", "lon"], SSTH_masked.data),
            "TH_850_masked": (["time", "lat", "lon"], TH_850_masked.data),
            "SHFLX_masked": (["time", "lat", "lon"], SHFLX_masked.data),
            "LHFLX_masked": (["time", "lat", "lon"], LHFLX_masked.data),
        },
        coords={
            "time": CAOI_masked.time,
            "lat": CAOI_masked.lat,
            "lon": CAOI_masked.lon,
        },
    )
    print('Saving data')

    # Save the dataset to a netCDF file
    output_path = f"{base_path}/{region}/"
    os.makedirs(output_path, exist_ok=True)
    file_path = f"{output_path}data_{timeframe}"

    # Check if file exists and if you have write permissions
    if os.path.exists(file_path):
        if os.access(file_path, os.W_OK):
            # Remove the file before overwriting (to avoid permission errors)
            os.remove(file_path)
        else:
            raise PermissionError(f"No write permission for {file_path}")
    ds_to_save.to_netcdf(file_path, mode='w', format='NETCDF4', engine='netcdf4')

def read_ensemble_data(region, timeframe, member_ids):
    ds_list = []
    for member_id in member_ids:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/data_{timeframe}"
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
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}"
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
    # Calculate SLP as hPa and add to dataset
    slp = all_members.PSL / 100  # Convert Pa to hPa
    slp.attrs["units"] = "hPa"
    all_members["SLP"] = slp
    # Calculate SLP anomalies
    slp_ano = all_members.SLP - all_members.SLP.mean(dim='time', skipna=True)
    slp_ano.attrs["units"] = "hPa"
    all_members["SLP_ano"] = slp_ano
    return all_members

def read_ensemble_data_ice(region, timeframe, member_ids):
    ds_list = []
    for member_id in member_ids:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/data_ice_{timeframe}"
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


def calculate_quantiles(region, timeframe):
    """Calculate 10th percentile of anomalies of DJF (in box) for quantiles"""
    CAOI_djf, CAOI_djf_ano = read_masked_CAOI_data(region, timeframe, member_ids)
    CAOI_djf_mean = calculate_space_mean(CAOI_djf)
    CAOI_djf_ano_mean = calculate_space_mean(CAOI_djf_ano)
    q10 = np.quantile(CAOI_djf_mean, 0.1)
    q90 = np.quantile(CAOI_djf_mean, 0.9)
    q10_ano = np.quantile(CAOI_djf_ano_mean, 0.1)
    q90_ano = np.quantile(CAOI_djf_ano_mean, 0.9)
    return q10, q90, q10_ano, q90_ano

def percentile_mask_time(region, timeframe, member_id):
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/clim_CAOI_{timeframe}_sorted_timesteplist_{region}"
    sorted_CAOI = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value"])
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/data_{timeframe}", decode_times=True)
    time_coord = ds.time

    if np.issubdtype(time_coord.dtype, np.datetime64):
        sorted_CAOI["time"] = pd.to_datetime(sorted_CAOI["time"], format="%Y-%m-%d-%H")
    else:
        import cftime
        sorted_CAOI["time"] = sorted_CAOI["time"].apply(
            lambda x: cftime.DatetimeNoLeap(*map(int, x.split("-")))
        )
    mask = time_coord.isin(sorted_CAOI["time"])
    return mask

def event_mask_time(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    """
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/data_{timeframe}", decode_times=True)
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

def event_mask_time_corrected(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    """
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_corrected_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/data_{timeframe}", decode_times=True)
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

def event_mask_time_precip(region, timeframe, member_id):
    """
    Create a mask for the CAOI_djf data based on timesteps of summary events for a specific ensemble member.
    """
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
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
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
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
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
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
    file_path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/summary_events_DJF90_{timeframe}_overstep_24_CAOI_{region}_filtered"
    summary_events = pd.read_csv(file_path, sep=" ", header=None, names=["time", "value", "event_ID", "duration"])

    # Open the member's dataset to get the time coordinate
    ds = xr.open_dataset(f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}/precip_data_time_extended{timeframe}", decode_times=True)
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
    # return xr.concat(masked_list, dim='member')


def get_percentile_in_box_ensemble(ds_ens, region, timeframe, varname, member_ids):
    """ Returns DataArray with dims (time, lat, lon) for all ensemble members in the box"""
    ### old version
    # def get_percentile_in_box(ds,region, timeframe, member_id):
    #     """Get the percentile timesteps in the box (lat, lon)"""
    #     lon_mask, lat_mask = box_mask(region)
    #     time_mask = percentile_mask_time(region, timeframe, member_id)
    #     ds_box = ds.sel(lon=lon_mask, lat=lat_mask, time =time_mask)
    #     return ds_box
    # def box_mask(region):
    #     lon_mask = slice(box[region][0], box[region][1])
    #     lat_mask = slice(box[region][2], box[region][3])
    #     return lon_mask, lat_mask
    # masked_list = []
    # for member_id in member_ids:
    #     ds_member = ds_ens.sel(member=member_id)
    #     ds_box = get_percentile_in_box(ds_member, region, timeframe, member_id)
    #     masked_list.append(ds_box)
    # return xr.concat(masked_list, dim=pd.Index(member_ids, name="member"))

    # ds_time_masked = get_masked_ensemble(ds_ens, region, timeframe, varname, percentile_mask_time, member_ids)
    # ds_box_masked = ds_time_masked.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3]))
    

    ds_box_masked = ds_ens.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3]))
    ds_time_masked = get_masked_ensemble(ds_box_masked, region, timeframe, varname, percentile_mask_time, member_ids)
    
    return ds_time_masked

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


# Histogram present vs future
def plot_histogram(region):

    CAOI_present = get_percentile_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    CAOI_future = get_percentile_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)

    CAOI_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    CAOI_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)

    CAOI_present_mean_space = calculate_space_mean(CAOI_present)
    CAOI_future_mean_space = calculate_space_mean(CAOI_future)

    CAOI_present_mean_space_time = np.nanmean(CAOI_present_mean_space)
    CAOI_future_mean_space_time = np.nanmean(CAOI_future_mean_space)

    q10_present = quantiles_dict[f'{region}_present']
    q10_future = quantiles_dict[f'{region}_future']
    print(f'present\nq10_ano_DJF = {q10_present:.2f}')
    print(f'future\nq10_ano_DJF = {q10_future:.2f}')


    plt.figure(figsize=(10, 5), layout='constrained')
    counts_present, bins_present, _ = plt.hist(
        CAOI_present_mean_space, bins=30, alpha=0.6, 
        label=f'present', 
        color='blue', density=True
    )
    plt.axvline(x=np.nanmean(CAOI_present_mean_space), color='darkblue', linestyle='dashed', linewidth=2, label = f'mean = {CAOI_present_mean_space_time:.2f}')
    
    counts_future, bins_future, _ = plt.hist(
        CAOI_future_mean_space, bins=30, alpha=0.6, 
        label=f'future', 
        color='orange', density=True
    )
    plt.axvline(x=np.nanmean(CAOI_future_mean_space), color='darkorange', linestyle='dashed', linewidth=2, label = f'mean = {CAOI_future_mean_space_time:.2f}')
    
    # plt.title(f"Histogram of CAOI for DJF - {region}")
    plt.xlim(2, 16)
    plt.ylim(0, 1)
    plt.xlabel("CAOI (K)", size=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel("frequency", size=20, labelpad=10)
    plt.legend(fontsize=20)
    plt.savefig(f"{fig_path}/CAOI_histogram_{region}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

    # Calculate the area under the bins
    area_present = np.sum(counts_present * np.diff(bins_present))
    area_future = np.sum(counts_future * np.diff(bins_future))
    print(f"Area under the Present histogram: {area_present}")
    print(f"Area under the Future histogram: {area_future}")

# Boxplot of DJF present vs future
def plot_boxplot(region):
    CAOI_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    CAOI_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)

    CAOI_present_mean_space = calculate_space_mean(CAOI_present)
    CAOI_future_mean_space = calculate_space_mean(CAOI_future)

    CAOI_present_mean_space_time = np.nanmean(CAOI_present_mean_space)
    CAOI_future_mean_space_time = np.nanmean(CAOI_future_mean_space)

    plt.figure(figsize=(6, 10))
    plt.boxplot([CAOI_present_mean_space, CAOI_future_mean_space], labels=[f'present\nmean = {CAOI_present_mean_space_time:.2f}', f'future\nmean = {CAOI_future_mean_space_time:.2f}'])
    # plt.title(f"Boxplot of CAOI for DJF - {region}")
    plt.ylim(2, 16)
    # print(plt.ylim())
    plt.ylabel("CAOI (K)", size=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=16)
    plt.savefig(f"{fig_path}/CAOI_boxplot_{region}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

def cmap_pcolormesh(cmap_name, boundaries, one_way=True):
    # Check if it's a cmcrameri colormap first
    if hasattr(cmc, cmap_name):
        base_cmap = getattr(cmc, cmap_name)
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

def plot_single_map_subplot(data_2d, lons, lats, region, timeframe_name, extent, cmap_name, boundaries, varname, 
                   labname=None, cbarticks=None, extend='neither', save_path=False, lag=0, 
                   significance_mask=None, ax=None, add_colorbar=True):
    """
    Plot a single map with optional external axis support for subplot creation.
    
    Parameters:
    -----------
    ax : matplotlib.axes or None
        If provided, plot on this axis instead of creating new figure
    add_colorbar : bool
        Whether to add colorbar (set to False when using as subplot)
    
    Returns:
    --------
    im : matplotlib plot object
        The plot object for colorbar creation
    """ 
    # Determine if this is a standalone plot or subplot
    standalone_plot = (ax is None)
    # set up map projection, event mask and extent depending on variable and region
    if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
        # For precipitation, evaporation, and P-E, use the map extent for the box
        if region == "NA":
            # For North Atlantic, use a different projection
            map_projection = ccrs.AzimuthalEquidistant(central_longitude=-55, central_latitude=50)
            # extent = [map_extent_precip[region][0]-360, map_extent_precip[region][1]-15, map_extent_precip[region][2], map_extent_precip[region][3]]
            extent = map_extent_precip[region]
        elif region == 'JS':
            map_projection = ccrs.AzimuthalEquidistant(central_longitude=160, central_latitude=50)
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
            SLP = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SLP_ano", event_mask_time_func, member_ids)
            SST = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
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
        else:
            # Normal case
            sea_ice = ice_masked_ensemble.sel(
                lon=slice(area_slice[0], area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
        print('Sea ice data shape:', sea_ice.shape)
        print(sea_ice.lat.max(), sea_ice.lat.min(), sea_ice.lon.max(), sea_ice.lon.min())
        

    elif timeframe_name == "10th_future" or timeframe_name == 'future':
        ds_name = 'future'
        # SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SLP = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SLP_ano", event_mask_time_func, member_ids)
            SST = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
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
        else:
            # Normal case
            sea_ice = ice_masked_ensemble.sel(
                lon=slice(area_slice[0], area_slice[1]), 
                lat=slice(area_slice[2], area_slice[3])
            )
 
        print('Sea ice data shape:', sea_ice.shape)
        print(sea_ice.lat.max(), sea_ice.lat.min(), sea_ice.lon.max(), sea_ice.lon.min())

    elif timeframe_name == "10th_diff" or timeframe_name == 'diff':
        # SST_present = get_masked_ensemble(data_dict[f"{region}_present"], region, "present", "SST", event_mask_time_func, member_ids)
        # SST_future = get_masked_ensemble(data_dict[f"{region}_future"], region, "future", "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
            SLP_present = get_masked_ensemble(data_dict_precip[f"{region}_present"], region, "present", "SLP_ano", event_mask_time_func, member_ids)
            SLP_future = get_masked_ensemble(data_dict_precip[f"{region}_future"], region, "future", "SLP_ano", event_mask_time_func, member_ids)
            SST_present = get_masked_ensemble(data_dict_precip[f"{region}_present"], region, "present", "SST", event_mask_time_func, member_ids)
            SST_present = SST_present.where(SST_present > 273)  # Mask out SST values below 0°C
            SST_future = get_masked_ensemble(data_dict_precip[f"{region}_future"], region, "future", "SST", event_mask_time_func, member_ids)
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
    # Create figure and axis only if not provided
    if standalone_plot:
        fig, ax = plt.subplots(1, 1, figsize=get_figsize_from_map_extent(extent, base_height=7), 
                              subplot_kw=dict(projection=map_projection), layout='constrained')
    else:
        fig = ax.figure
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.coastlines(color='black')
    ax.set_facecolor("white")

    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 20, "color": "k"}
    gl.ylabel_style = {"size": 20, "color": "k"}
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
    
    if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
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
            SST_mean = SST_mean.where(mask)
            SST_contour = ax.contour(lons_plot, lats, SST_mean, levels = np.arange(275, 305, 5), colors = 'black', linewidths = 2, transform=ccrs.PlateCarree())
            ax.clabel(SST_contour, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
            ax.contour(lons_plot, lats, calculate_time_mean(SLP), levels=np.arange(-6, 0, 2), colors='dimgrey', linestyles = 'dashed',linewidths=2, transform=ccrs.PlateCarree(), alpha=1)
            ax.contour(lons_plot, lats, calculate_time_mean(SLP), levels=np.arange(0, 8, 2), colors='dimgrey', linewidths=2, transform=ccrs.PlateCarree(), alpha=1)
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
            SST_contour = ax.contour(
                SST_mean.lon, SST_mean.lat, SST_mean,
                levels=np.arange(275, 310, 5),
                colors='black', linewidths=2,
                transform=ccrs.PlateCarree()
            )
            ax.clabel(SST_contour, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
        
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
            SST_contour_neg = ax.contour(lons_plot, lats, calculate_time_mean(SST_future)-calculate_time_mean(SST_present), levels = [-4, -3, -2, -1], colors = 'black', negative_linestyles = 'dashed', linewidths = 2, transform=ccrs.PlateCarree())
            SST_contour_pos = ax.contour(lons_plot, lats, calculate_time_mean(SST_future)-calculate_time_mean(SST_present), levels = [1, 2, 3, 4, 5, 6], colors = 'black', negative_linestyles = 'dashed', linewidths = 2, transform=ccrs.PlateCarree())
            ax.clabel(SST_contour_neg, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
            ax.clabel(SST_contour_pos, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
        if varname == 'CAOI':
            # Get SST masked with the same event mask and extent as CAOI present and future
            SST_present_event_masked = get_masked_ensemble(
                data_dict[f"{region}_present"], region, "present", "SST", event_mask_time, member_ids
            )
            SST_present_event_masked = SST_present_event_masked.where(SST_present_event_masked > 273)
            SST_future_event_masked = get_masked_ensemble(
                data_dict[f"{region}_future"], region, "future", "SST", event_mask_time, member_ids
            )
            SST_future_event_masked = SST_future_event_masked.where(SST_future_event_masked > 273)
            # Select the same extent as CAOI
            SST_present_event_masked = SST_present_event_masked.sel(
                lon=slice(extent[0], extent[1]), lat=slice(extent[2], extent[3])
            )
            SST_future_event_masked = SST_future_event_masked.sel(
                lon=slice(extent[0], extent[1]), lat=slice(extent[2], extent[3])
            )
            # Calculate the time mean for contour plotting
            SST_present_mean = calculate_time_mean(SST_present_event_masked)
            SST_future_mean = calculate_time_mean(SST_future_event_masked)
            # Plot contour lines (e.g., every 2K from 273 to 305)
            SST_contour_neg = ax.contour(
                SST_present_mean.lon, SST_present_mean.lat, 
                SST_future_mean - SST_present_mean,
                levels=[-4, -3, -2, -1],
                colors='black', linestyles='dashed', linewidths=2,
                transform=ccrs.PlateCarree()
            )
            SST_contour_pos = ax.contour(
                SST_present_mean.lon, SST_present_mean.lat, 
                SST_future_mean - SST_present_mean,
                levels=[1, 2, 3, 4, 5, 6],
                colors='black', linewidths=2,
                transform=ccrs.PlateCarree()
            )
            ax.clabel(SST_contour_neg, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
            ax.clabel(SST_contour_pos, inline=True, fontsize=16, fmt='%d', colors='black', use_clabeltext=True)
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
    if significance_mask is not None and ('diff' in timeframe_name):
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
    # Only add colorbar and save for standalone plots
    if standalone_plot and add_colorbar:
        divider = make_axes_locatable(ax)
        shrink = get_colorbar_shrink(extent)
        if extent == map_extent[region]:
            shrink = 0.9
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal', shrink=shrink, 
                           pad=0.05, extend=extend)
        cbar.set_label(labname if labname else varname, size=20)
        cbar.ax.tick_params(labelsize=16)
        if cbarticks is not None:
            cbar.set_ticks(cbarticks)
        else:
            cbar.set_ticks(boundaries)

    if lag != 0:
        lag_name = f'_plus{lag}'
    else:
        lag_name = ''

    if standalone_plot:
        if save_path == True:
            save_path = f"{fig_path}/{varname}_{timeframe_name}_{region}{lag_name}.png"
            plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
            print(f"Saved plot to: {save_path}")

            if varname not in ['CAOI_8Kfrequency_10th']:
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
                print('Mean over box:', calculate_space_mean(data_2d_box))

        plt.show()
        plt.close()

    return im


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
            SLP = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SLP_ano", event_mask_time_func, member_ids)
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
            SLP = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SLP_ano", event_mask_time_func, member_ids)
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
            SLP_present = get_masked_ensemble(data_dict_precip[f"{region}_present"], region, "present", "SLP_ano", event_mask_time_func, member_ids)
            SLP_future = get_masked_ensemble(data_dict_precip[f"{region}_future"], region, "future", "SLP_ano", event_mask_time_func, member_ids)
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
    if significance_mask is not None and ('diff' in timeframe_name):
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

def plot_map_djf(region):
    CAOI_djf_present = data_dict[f'{region}_present'].CAOI.mean(dim='member', skipna=True)
    CAOI_djf_future = data_dict[f'{region}_future'].CAOI.mean(dim='member', skipna=True)

    CAOI_djf_present_mean_time = calculate_time_mean(CAOI_djf_present)
    CAOI_djf_future_mean_time = calculate_time_mean(CAOI_djf_future)

    CAOI_djf_diff_mean_time = CAOI_djf_future_mean_time - CAOI_djf_present_mean_time

    lonmin, lonmax = CAOI_djf_present.lon.min(), CAOI_djf_present.lon.max()
    latmin, latmax = CAOI_djf_present.lat.min(), CAOI_djf_present.lat.max()
    box = [lonmin, lonmax, latmin, latmax]

    lons = CAOI_djf_present.lon
    lats = CAOI_djf_present.lat

    plot_single_map(
        data_2d=CAOI_djf_present_mean_time,
        lons =lons,
        lats =lats,
        region=region,
        timeframe_name='djf_present',
        extent= box,
        cmap_name='plasma',
        boundaries=np.arange(0, 13, 1),
        varname='CAOI',
        labname='CAOI (K)',
        extend='max',
        save_path=True
    )

    plot_single_map(
        data_2d=CAOI_djf_future_mean_time,
        lons =lons,
        lats =lats,
        region=region,
        timeframe_name='djf_future',
        extent= box,
        cmap_name='plasma',
        boundaries=np.arange(0, 13, 1),
        varname='CAOI',
        labname='CAOI (K)',
        extend='max',
        save_path=True
    )

    plot_single_map(
        data_2d=CAOI_djf_diff_mean_time,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='djf_diff',
        extent= box,
        cmap_name='RdBu_r',
        boundaries=np.arange(-3, 3.5, 0.5),
        cbarticks=np.arange(-3, 4, 1),
        varname='CAOI',
        labname='Δ CAOI (K)',
        extend='both',
        save_path=True
    )

def plot_map_10th(region, member_ids=member_ids, varname='CAOI', lag = 0):
    lag = lag
    # Get the configuration for this variable
    if varname not in VAR_CONFIG:
        raise ValueError(f"Variable '{varname}' not found in VAR_CONFIG. Available variables: {list(VAR_CONFIG.keys())}")
    
    if varname in ['PRECIP', 'EVAP', 'EP', 'SST']:
        # For precipitation, evaporation, and P-E, use the read_ensemble_data_precip function
        data_dict_v = data_dict_precip
        if lag == 0:
            event_mask_time_func = event_plus0_mask_time
        elif lag == 1:
            event_mask_time_func = event_plus1_mask_time
        elif lag == 2:
            event_mask_time_func = event_plus2_mask_time
        
    else:
        data_dict_v = data_dict
        event_mask_time_func = event_mask_time
        

    config = VAR_CONFIG[varname]
    diff_config = VAR_CONFIG_DIFF[varname]

    

    CAOI_10th_present = get_masked_ensemble(data_dict_v[f'{region}_present'], region, "present", varname, event_mask_time_func, member_ids)
    CAOI_10th_future = get_masked_ensemble(data_dict_v[f'{region}_future'], region, "future", varname, event_mask_time_func, member_ids)

    # if lag == 0:
    #     CAOI_10th_present = get_masked_ensemble(data_dict_v[f'{region}_present'], region, "present", varname, percentile_mask_time, member_ids)
    #     CAOI_10th_future = get_masked_ensemble(data_dict_v[f'{region}_future'], region, "future", varname, percentile_mask_time, member_ids)
    # elif lag == 1:
    #     CAOI_10th_present = get_masked_ensemble(data_dict_v[f'{region}_present'], region, "present", varname, event_plus1_mask_time, member_ids)
    #     CAOI_10th_future = get_masked_ensemble(data_dict_v[f'{region}_future'], region, "future", varname, event_plus1_mask_time, member_ids)
    # elif lag == 2:
    #     CAOI_10th_present = get_masked_ensemble(data_dict_v[f'{region}_present'], region, "present", varname, event_plus2_mask_time, member_ids)
    #     CAOI_10th_future = get_masked_ensemble(data_dict_v[f'{region}_future'], region, "future", varname, event_plus2_mask_time, member_ids)
    # Now CAOI_10th_present has dims (member, time, lat, lon)
    CAOI_10th_present_mean_time = calculate_time_mean(CAOI_10th_present)
    CAOI_10th_future_mean_time = calculate_time_mean(CAOI_10th_future)
    CAOI_10th_diff_mean_time = CAOI_10th_future_mean_time - CAOI_10th_present_mean_time

    # Use your plotting function as before
    lonmin, lonmax = CAOI_10th_present.lon.min(), CAOI_10th_present.lon.max()
    latmin, latmax = CAOI_10th_present.lat.min(), CAOI_10th_present.lat.max()
    box = [lonmin, lonmax, latmin, latmax]

    lons = CAOI_10th_present.lon
    lats = CAOI_10th_present.lat

    # if varname == 'PRECIP' or varname == 'EVAP' or varname == 'EP':
    #     box = map_extent_precip[region]

    # box = map_extent_precip[region]

    plot_single_map(
        data_2d=CAOI_10th_present_mean_time,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='10th_present',
        extent=box,
        cmap_name=config['cmap_name'],
        boundaries=config['boundaries'],
        cbarticks=config['cbarticks'],
        varname=varname,
        labname=config['labname'],
        extend='max',
        save_path=True,
        lag=lag
    )

    plot_single_map(
        data_2d=CAOI_10th_future_mean_time,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='10th_future',
        extent=box,
        cmap_name=config['cmap_name'],
        boundaries=config['boundaries'],
        cbarticks=config['cbarticks'],
        varname=varname,
        labname=config['labname'],
        extend='max',
        save_path=True,
        lag=lag
    )

    plot_single_map(
        data_2d=CAOI_10th_diff_mean_time,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='10th_diff',
        extent=box,
        cmap_name=diff_config['cmap_name'],
        boundaries=diff_config['boundaries'],
        cbarticks=diff_config['cbarticks'],
        varname=varname,
        labname=diff_config['labname'],
        extend='both',
        save_path=True,
        lag=lag
    )

def compute_lag1_autocorr(x):
    """Compute lag-1 autocorrelation ignoring NaNs, safely."""
    x = x - np.nanmean(x)
    num = np.nansum(x[:-1] * x[1:])
    denom = np.nansum(x**2)
    if denom == 0:
        return np.nan  # or 0 if you want to assume no autocorrelation
    return num / denom


def compute_ztest_with_neff(CAOI_present, CAOI_future):
    time, lat, lon = CAOI_present.shape

    x1 = np.nansum(CAOI_present > 8, axis=0)
    x2 = np.nansum(CAOI_future > 8, axis=0)

    n1 = np.sum(~np.isnan(CAOI_present), axis=0)
    n2 = np.sum(~np.isnan(CAOI_future), axis=0)

    p1 = x1 / n1
    p2 = x2 / n2
    p_pool = (x1 + x2) / (n1 + n2)

    # Compute lag-1 autocorrelation and neff
    rho1 = np.full((lat, lon), np.nan)
    rho2 = np.full((lat, lon), np.nan)

    for i in range(lat):
        for j in range(lon):
            ts1 = (CAOI_present[:, i, j] > 8).astype(float)
            ts2 = (CAOI_future[:, i, j] > 8).astype(float)

            ts1[~np.isfinite(ts1)] = np.nan
            ts2[~np.isfinite(ts2)] = np.nan

            if np.sum(np.isfinite(ts1)) > 2:
                rho1[i, j] = compute_lag1_autocorr(ts1)
            else:
                rho1[i, j] = 0  # fallback: assume no autocorrelation

            if np.sum(np.isfinite(ts2)) > 2:
                rho2[i, j] = compute_lag1_autocorr(ts2)
            else:
                rho2[i, j] = 0

    # Bound rho to avoid instability
    rho1 = np.clip(rho1, -0.99, 0.99)
    rho2 = np.clip(rho2, -0.99, 0.99)

    # Effective sample size
    neff1 = n1 * (1 - rho1) / (1 + rho1)
    neff2 = n2 * (1 - rho2) / (1 + rho2)

    # Standard error and Z-test
    with np.errstate(divide='ignore', invalid='ignore'):
        se = np.sqrt(p_pool * (1 - p_pool) * (1 / neff1 + 1 / neff2))
        z = (p2 - p1) / se
        p_values = 2 * (1 - norm.cdf(np.abs(z)))

    significance_mask = p_values < 0.05

    return p1 * 100, p2 * 100, (p2 - p1) * 100, significance_mask


def plot_map_frequency(region):

    CAOI_djf_present = data_dict[f'{region}_present'].CAOI
    CAOI_djf_future = data_dict[f'{region}_future'].CAOI

    lonmin, lonmax = CAOI_djf_present.lon.min(), CAOI_djf_present.lon.max()
    latmin, latmax = CAOI_djf_present.lat.min(), CAOI_djf_present.lat.max()
    box = [lonmin, lonmax, latmin, latmax]
    lons = CAOI_djf_present.lon
    lats = CAOI_djf_present.lat

    # Frequency of CAOI > 8 K
    # Count for each grid point how many times the CAOI is above 8 K, ignoring NaNs
    CAOI_present_count = np.nansum(CAOI_djf_present > 8, axis=(0, 1))  # Sum over time and member dimensions
    CAOI_future_count = np.nansum(CAOI_djf_future > 8, axis=(0, 1))  # Sum over time and member dimensions

    # Divide by the total number of valid timesteps (non-NaN) to get the frequency
    CAOI_present_freq = CAOI_present_count / np.sum(~np.isnan(CAOI_djf_present), axis=(0, 1)) * 100
    print('max present', np.max(CAOI_present_freq))
    print('mean present', np.nanmean(CAOI_present_freq))
    CAOI_future_freq = CAOI_future_count / np.sum(~np.isnan(CAOI_djf_future), axis=(0, 1)) * 100
    print('max future', np.max(CAOI_future_freq))
    print('mean future', np.nanmean(CAOI_future_freq))
    CAOI_diff_freq = CAOI_future_freq - CAOI_present_freq
    print('max diff', np.max(CAOI_diff_freq))
    print('mean diff', np.nanmean(CAOI_diff_freq))

    # # Divide by the total number of valid timesteps (non-NaN) to get the frequency
    # CAOI_present_freq = CAOI_present_count / np.sum(~np.isnan(CAOI_present), axis=(0)) * 100
    # print('max present', np.max(CAOI_present_freq))
    # print('mean present', np.nanmean(CAOI_present_freq))
    # CAOI_future_freq = CAOI_future_count / np.sum(~np.isnan(CAOI_future), axis=(0)) * 100
    # print('max future', np.max(CAOI_future_freq))
    # print('mean future', np.nanmean(CAOI_future_freq))
    # CAOI_diff_freq = CAOI_future_freq - CAOI_present_freq
    # print('max diff', np.max(CAOI_diff_freq))
    # print('mean diff', np.nanmean(CAOI_diff_freq))

    # Divide by the total number of timesteps to get the frequency
    CAOI_present_freq = CAOI_present_count / (len(CAOI_djf_present.time)*len(CAOI_djf_present.member)) * 100
    CAOI_future_freq = CAOI_future_count / (len(CAOI_djf_future.time)*len(CAOI_djf_future.member)) * 100
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

    plot_single_map(
        data_2d=CAOI_present_freq,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='present',
        extent=box,
        cmap_name='Oranges',
        boundaries=np.arange(0, 55, 5),
        varname='CAOI_8Kfrequency_djf',
        labname='frequency exceedance 8K (%)',
        extend='max',
        save_path=True
    )

    plot_single_map(
        data_2d=CAOI_future_freq,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='future',
        extent=box,
        cmap_name='Oranges',
        boundaries=np.arange(0, 55, 5),
        varname='CAOI_8Kfrequency_djf',
        labname='frequency exceedance 8K (%)',
        extend='max',
        save_path=True
    )

    plot_single_map(
        data_2d=CAOI_diff_freq,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='diff',
        extent=box,
        cmap_name='RdBu_r',
        boundaries=np.arange(-20, 22.5, 2.5),
        cbarticks=np.arange(-20, 22.5, 5),
        varname='CAOI_8Kfrequency_djf',
        labname='Δ frequency exceedance 8K (%)',
        extend='both',
        save_path=True
    )

def plot_map_frequency_10th(region):

    CAOI_present = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', event_mask_time, member_ids)
    CAOI_future = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', event_mask_time, member_ids)

    lonmin, lonmax = CAOI_present.lon.min(), CAOI_present.lon.max()
    latmin, latmax = CAOI_present.lat.min(), CAOI_present.lat.max()
    box = [lonmin, lonmax, latmin, latmax]
    lons = CAOI_present.lon
    lats = CAOI_present.lat

    # Frequency of CAOI > 8 K
    # Count for each grid point how many times the CAOI is above 8 K, ignoring NaNs
    CAOI_present_count = np.nansum(CAOI_present > 8, axis=(0))  # Sum over time and member dimensions
    CAOI_future_count = np.nansum(CAOI_future > 8, axis=(0))  # Sum over time and member dimensions

    # # Divide by the total number of valid timesteps (non-NaN) to get the frequency
    # CAOI_present_freq = CAOI_present_count / np.sum(~np.isnan(CAOI_present), axis=(0)) * 100
    # print('max present', np.max(CAOI_present_freq))
    # print('mean present', np.nanmean(CAOI_present_freq))
    # CAOI_future_freq = CAOI_future_count / np.sum(~np.isnan(CAOI_future), axis=(0)) * 100
    # print('max future', np.max(CAOI_future_freq))
    # print('mean future', np.nanmean(CAOI_future_freq))
    # CAOI_diff_freq = CAOI_future_freq - CAOI_present_freq
    # print('max diff', np.max(CAOI_diff_freq))
    # print('mean diff', np.nanmean(CAOI_diff_freq))

    # Divide by the total number of timesteps to get the frequency
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

    plot_single_map(
        data_2d=CAOI_present_freq,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='present',
        extent=box,
        cmap_name='bilbao',
        boundaries=np.arange(0, 110, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='frequency exceedance 8K [%]',
        extend='neither',
        save_path=True
    )

    plot_single_map(
        data_2d=CAOI_future_freq,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='future',
        extent=box,
        cmap_name='bilbao',
        boundaries=np.arange(0, 110, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='frequency exceedance 8K [%]',
        extend='neither',
        save_path=True
    )

    plot_single_map(
        data_2d=CAOI_diff_freq,
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='diff',
        extent=box,
        cmap_name='vik',
        boundaries=np.arange(-30, 35, 5),
        cbarticks=np.arange(-30, 40, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='Δ frequency exceedance 8K [%]',
        extend='both',
        save_path=True
    )

def plot_map_fluxes_10th(region):
    SHFLX_present = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'SHFLX_masked', event_mask_time, member_ids)
    SHFLX_future = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'SHFLX_masked', event_mask_time, member_ids)
    LHFLX_present = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'LHFLX_masked', event_mask_time, member_ids)
    LHFLX_future = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'LHFLX_masked', event_mask_time, member_ids)

    box = [data_dict[f'{region}_present'].lon.min(), data_dict[f'{region}_present'].lon.max(), data_dict[f'{region}_present'].lat.min(), data_dict[f'{region}_present'].lat.max()]
    lons = data_dict[f'{region}_present'].lon
    lats = data_dict[f'{region}_present'].lat

    dict_boundaries = {
        'SHFLX_JS': np.arange(0, 361, 20),
        'SHFLX_NA': np.arange(0, 481, 24),
        'LHFLX_JS': np.arange(0, 561, 40),
        'LHFLX_NA': np.arange(0, 801, 40),
        'THFLX_JS': np.arange(0, 801, 50),
        'THFLX_NA': np.arange(0, 1201, 60),
        'SHFLX_diff': np.arange(-150, 151, 15),
        'LHFLX_diff': np.arange(-100, 101, 10),
        'THFLX_diff': np.arange(-200, 201, 20)
    }
    dict_cbarticks = {
        'SHFLX_JS': np.arange(0, 361, 40),
        'SHFLX_NA': np.arange(0, 481, 48),
        'LHFLX_JS': np.arange(0, 561, 80),
        'LHFLX_NA': np.arange(0, 801, 80),
        'THFLX_JS': np.arange(0, 801, 100),
        'THFLX_NA': np.arange(0, 1201, 120),
        'SHFLX_diff': np.arange(-150, 151, 30),
        'LHFLX_diff': np.arange(-100, 101, 20),
        'THFLX_diff': np.arange(-200, 201, 40)
    }


    plot_single_map(
        data_2d=calculate_time_mean(SHFLX_present),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='present',
        extent=box,
        cmap_name='Reds',
        boundaries=dict_boundaries[f'SHFLX_NA'],
        cbarticks=dict_cbarticks[f'SHFLX_NA'],
        varname='FLX_SH',
        labname='surface sensible heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d=calculate_time_mean(SHFLX_future),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='future',
        extent=box,
        cmap_name='Reds',
        boundaries=dict_boundaries[f'SHFLX_NA'],
        cbarticks=dict_cbarticks[f'SHFLX_NA'],
        varname='FLX_SH',
        labname='surface sensible heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d = calculate_time_mean(SHFLX_future)-calculate_time_mean(SHFLX_present),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='diff',
        extent=box,
        cmap_name='RdBu_r',
        boundaries=dict_boundaries[f'SHFLX_diff'],
        cbarticks=dict_cbarticks[f'SHFLX_diff'],
        varname='FLX_SH',
        labname='Δ surface sensible heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d = calculate_time_mean(LHFLX_present),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='present',
        extent=box,
        cmap_name='Purples',
        boundaries=dict_boundaries[f'LHFLX_NA'],
        cbarticks=dict_cbarticks[f'LHFLX_NA'],
        varname='FLX_LH',
        labname='surface latent heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d = calculate_time_mean(LHFLX_future),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='future',
        extent=box,
        cmap_name='Purples',
        boundaries=dict_boundaries[f'LHFLX_NA'],
        cbarticks=dict_cbarticks[f'LHFLX_NA'],
        varname='FLX_LH',
        labname='surface latent heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d = calculate_time_mean(LHFLX_future)-calculate_time_mean(LHFLX_present),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='diff',
        extent=box,
        cmap_name='RdBu_r',
        boundaries=dict_boundaries[f'LHFLX_diff'],
        cbarticks=dict_cbarticks[f'LHFLX_diff'],
        varname='FLX_LH',
        labname='Δ surface latent heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d=calculate_time_mean(SHFLX_present+LHFLX_present),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='present',
        extent=box,
        cmap_name='Reds',
        boundaries=dict_boundaries[f'THFLX_NA'],
        cbarticks=dict_cbarticks[f'THFLX_NA'],
        varname='FLX_TH',
        labname='surface tubulent heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d=calculate_time_mean(SHFLX_future+LHFLX_future),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='future',
        extent=box,
        cmap_name='Reds',
        boundaries=dict_boundaries[f'THFLX_NA'],
        cbarticks=dict_cbarticks[f'THFLX_NA'],
        varname='FLX_TH',
        labname='surface tubulent heat flux (W/m2)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d = calculate_time_mean(SHFLX_future+LHFLX_future)-calculate_time_mean(SHFLX_present+LHFLX_present),
        lons=lons,
        lats=lats,
        region=region,
        timeframe_name='diff',
        extent=box,
        cmap_name='RdBu_r',
        boundaries=dict_boundaries[f'THFLX_diff'],
        cbarticks=dict_cbarticks[f'THFLX_diff'],
        varname='FLX_TH',
        labname='Δ surface tubulent heat flux (W/m2)',
        extend='both',
        save_path=True
    )

def plot_map_fluxes_contribution(region):
    for timeframe in ["present", "future"]:
        # Get DJF fluxes (all DJF times)
        SHFLX_djf = data_dict[f'{region}_{timeframe}'].SHFLX_masked
        LHFLX_djf = data_dict[f'{region}_{timeframe}'].LHFLX_masked

        # Get 10th percentile fluxes (masked times)
        SHFLX_10th = get_masked_ensemble(data_dict[f'{region}_{timeframe}'], region, timeframe, 'SHFLX_masked', event_mask_time, member_ids)
        LHFLX_10th = get_masked_ensemble(data_dict[f'{region}_{timeframe}'], region, timeframe, 'LHFLX_masked', event_mask_time, member_ids)

        # For each gridpoint: sum over time for 10th and DJF, then ratio
        SHFLX_10th_sum = SHFLX_10th.sum(dim=['time'], skipna=True)
        SHFLX_djf_sum = SHFLX_djf.sum(dim=['time'], skipna=True)

        # For LHFLX: only sum positive values
        LHFLX_10th_positive = LHFLX_10th.where(LHFLX_10th > 0)
        LHFLX_djf_positive = LHFLX_djf.where(LHFLX_djf > 0)
        LHFLX_10th_sum = LHFLX_10th_positive.sum(dim=['time'], skipna=True)
        LHFLX_djf_sum = LHFLX_djf_positive.sum(dim=['time'], skipna=True)

        # Print number of positive values for each gridpoint
        LHFLX_10th_positive_count = LHFLX_10th_positive.count(dim=['time'])
        LHFLX_10th_count = LHFLX_10th.count(dim=['time'])
        LHFLX_djf_positive_count = LHFLX_djf_positive.count(dim=['time'])
        print(f"LHFLX 10th positive count ({region}, {timeframe}):\n", LHFLX_10th_positive_count.values)
        print(f"LHFLX DJF positive count ({region}, {timeframe}):\n", LHFLX_djf_positive_count.values)

        # Contribution as percentage (per gridpoint)
        SHFLX_contrib_map = (SHFLX_10th_sum / SHFLX_djf_sum * 100)
        LHFLX_contrib_map = (LHFLX_10th_sum / LHFLX_djf_sum * 100)

        # Plot SHFLX contribution map
        plot_single_map(
            data_2d=SHFLX_contrib_map,
            lons=SHFLX_contrib_map.lon,
            lats=SHFLX_contrib_map.lat,
            region=region,
            timeframe_name=f'flux_contrib_SHFLX_{timeframe}',
            extent=[SHFLX_contrib_map.lon.min(), SHFLX_contrib_map.lon.max(), SHFLX_contrib_map.lat.min(), SHFLX_contrib_map.lat.max()],
            cmap_name='Reds',
            boundaries=np.arange(0, 101, 10),
            varname='SHFLX_contribution',
            labname=f'10th percentile SHFLX / DJF SHFLX (%) ({timeframe})',
            extend='neither',
            save_path=True
        )

        # plot ratio of positive LHFLX counts
        LHFLX_contrib_count_map = (LHFLX_10th_positive_count / LHFLX_10th_count * 100)

        plot_single_map(
            data_2d=LHFLX_contrib_count_map,
            lons=LHFLX_contrib_map.lon,
            lats=LHFLX_contrib_map.lat,
            region=region,
            timeframe_name=f'flux_contrib_LHFLX_{timeframe}',
            extent=[LHFLX_contrib_map.lon.min(), LHFLX_contrib_map.lon.max(), LHFLX_contrib_map.lat.min(), LHFLX_contrib_map.lat.max()],
            cmap_name='plasma',
            boundaries=np.arange(99, 102, 1),
            varname='LHFLX_contribution',
            labname=f'10th percentile LHFLX / DJF LHFLX (%) ({timeframe})',
            extend='neither',
            save_path=False
        )
        # Plot LHFLX contribution map
        plot_single_map(
            data_2d=LHFLX_contrib_map,
            lons=LHFLX_contrib_map.lon,
            lats=LHFLX_contrib_map.lat,
            region=region,
            timeframe_name=f'flux_contrib_LHFLX_{timeframe}',
            extent=[LHFLX_contrib_map.lon.min(), LHFLX_contrib_map.lon.max(), LHFLX_contrib_map.lat.min(), LHFLX_contrib_map.lat.max()],
            cmap_name='Purples',
            boundaries=np.arange(0, 101, 10),
            varname='LHFLX_contribution',
            labname=f'10th percentile LHFLX / DJF LHFLX (%) ({timeframe})',
            extend='neither',
            save_path=True
        )

    # Compute and plot the change in contribution (future - present)
    # Get present and future contributions
    SHFLX_djf_present = data_dict[f'{region}_present'].SHFLX_masked
    SHFLX_djf_future = data_dict[f'{region}_future'].SHFLX_masked
    LHFLX_djf_present = data_dict[f'{region}_present'].LHFLX_masked
    LHFLX_djf_future = data_dict[f'{region}_future'].LHFLX_masked

    SHFLX_10th_present = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'SHFLX_masked', event_mask_time, member_ids)
    SHFLX_10th_future = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'SHFLX_masked', event_mask_time, member_ids)
    LHFLX_10th_present = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'LHFLX_masked', event_mask_time, member_ids)
    LHFLX_10th_future = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'LHFLX_masked', event_mask_time, member_ids)


    SHFLX_10th_sum_present = SHFLX_10th_present.sum(dim=['time'], skipna=True)
    SHFLX_10th_sum_future = SHFLX_10th_future.sum(dim=['time'], skipna=True)
    SHFLX_djf_sum_present = SHFLX_djf_present.sum(dim=['time'], skipna=True)
    SHFLX_djf_sum_future = SHFLX_djf_future.sum(dim=['time'], skipna=True)

    LHFLX_10th_sum_present = LHFLX_10th_present.sum(dim=['time'], skipna=True)
    LHFLX_10th_sum_future = LHFLX_10th_future.sum(dim=['time'], skipna=True)
    LHFLX_djf_sum_present = LHFLX_djf_present.sum(dim=['time'], skipna=True)
    LHFLX_djf_sum_future = LHFLX_djf_future.sum(dim=['time'], skipna=True)

    SHFLX_contrib_present = (SHFLX_10th_sum_present / SHFLX_djf_sum_present * 100)
    SHFLX_contrib_future = (SHFLX_10th_sum_future / SHFLX_djf_sum_future * 100)
    LHFLX_contrib_present = (LHFLX_10th_sum_present / LHFLX_djf_sum_present * 100)
    LHFLX_contrib_future = (LHFLX_10th_sum_future / LHFLX_djf_sum_future * 100)

    SHFLX_contrib_diff = SHFLX_contrib_future - SHFLX_contrib_present
    LHFLX_contrib_diff = LHFLX_contrib_future - LHFLX_contrib_present

    # Plot difference maps
    plot_single_map(
        data_2d=SHFLX_contrib_diff,
        lons=SHFLX_contrib_diff.lon,
        lats=SHFLX_contrib_diff.lat,
        region=region,
        timeframe_name='flux_contrib_SHFLX_diff',
        extent=[SHFLX_contrib_diff.lon.min(), SHFLX_contrib_diff.lon.max(), SHFLX_contrib_diff.lat.min(), SHFLX_contrib_diff.lat.max()],
        cmap_name='RdBu_r',
        boundaries=np.arange(-30, 35, 5),
        varname='SHFLX_contribution_diff',
        labname='Δ 10th percentile SHFLX / DJF SHFLX (%) (future - present)',
        extend='both',
        save_path=True
    )

    plot_single_map(
        data_2d=LHFLX_contrib_diff,
        lons=LHFLX_contrib_diff.lon,
        lats=LHFLX_contrib_diff.lat,
        region=region,
        timeframe_name='flux_contrib_LHFLX_diff',
        extent=[LHFLX_contrib_diff.lon.min(), LHFLX_contrib_diff.lon.max(), LHFLX_contrib_diff.lat.min(), LHFLX_contrib_diff.lat.max()],
        cmap_name='RdBu_r',
        boundaries=np.arange(-30, 35, 5),
        varname='LHFLX_contribution_diff',
        labname='Δ 10th percentile LHFLX / DJF LHFLX (%) (future - present)',
        extend='both',
        save_path=True
    )

def plot_caoi_contribution_bar(region):
    # PRESENT
    # -----------------------
    present_CAOI = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    present_SSTH = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    present_TH_850 = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)


    weights = np.cos(np.deg2rad(present_CAOI.lat))
    weights.name = "weights"
    present_CAOI = present_CAOI.weighted(weights)
    present_SSTH = present_SSTH.weighted(weights)
    present_TH_850 = present_TH_850.weighted(weights)




    # FUTURE
    # -----------------------
    future_CAOI = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    future_SSTH = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    future_TH_850 = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)

    weights = np.cos(np.deg2rad(future_CAOI.lat))
    weights.name = "weights"
    future_CAOI = future_CAOI.weighted(weights)
    future_SSTH = future_SSTH.weighted(weights)
    future_TH_850 = future_TH_850.weighted(weights)
    
    # Calculate means
    present_vals = [present_CAOI.mean(dim = ['time', 'lat', 'lon'], skipna=True), 
                   present_SSTH.mean(dim = ['time', 'lat', 'lon'], skipna=True),
                   present_TH_850.mean(dim = ['time', 'lat', 'lon'], skipna=True)]
    future_vals = [future_CAOI.mean(dim = ['time', 'lat', 'lon'], skipna=True),
                     future_SSTH.mean(dim = ['time', 'lat', 'lon'], skipna=True),
                     future_TH_850.mean(dim = ['time', 'lat', 'lon'], skipna=True)]
    # present_vals = list(map(np.nanmean, [present_CAOI_filtered, present_SSTH_filtered, present_TH_850_filtered]))
    # future_vals = list(map(np.nanmean, [future_CAOI_filtered, future_SSTH_filtered, future_TH_850_filtered]))
    print(present_vals)
    print(future_vals)


    x = [0, 1]  # Positions for SSTH and TH850 only
    width = 0.3
    labels = ['SST', '850', 'CAOI']
    if region == 'JS':
        ylims = (273, 288)
    elif region == 'NA':
        ylims = (280, 295)

    fig, ax1 = plt.subplots(figsize=(12, 9))

    # Plot SSTH and TH850
    ax1.bar([val - width/2 for val in x], present_vals[1:3], width=width,
            color=['dodgerblue', 'tomato'], label='present', alpha=0.7)
    ax1.bar([val + width/2 for val in x], future_vals[1:3], width=width,
            color=['skyblue', 'salmon'], label='future', alpha=0.7)

    ax1.set_ylabel(r"mean $\Theta$ (K)", fontsize=20)
    ax1.set_xticks([0, 1, 2])
    ax1.set_ylim(ylims)
    ax1.tick_params(axis='both', labelsize=16)
    ax1.grid(axis='y', linestyle='--', alpha=0.5)
    ax1.set_xticklabels(labels, fontsize=20)

    # Plot CAOI on secondary y-axis
    ax2 = ax1.twinx()
    ax2.bar(2 - width/2, present_vals[0], width=width,
            color='darkgreen', label='present CAOI', alpha=0.7)
    ax2.bar(2 + width/2, future_vals[0], width=width,
            color='lightgreen', label='future CAOI', alpha=0.7)
    ax2.set_ylabel("mean CAOI (K)", fontsize=20, color='darkgreen')
    ax2.set_ylim(0, 15)
    ax2.tick_params(axis='y', labelsize=16, colors='darkgreen')

    # Legends with fontsize
    ax1.legend(loc='upper left', fontsize=20)
    ax2.legend(loc='upper right', fontsize=20)

    # Save and show
    plt.savefig(f"{fig_path}/CAOI_contribution_bar_{region}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

def plot_caoi_contribution(region):
    if region == "JS":
        ylims_temp = (245, 305)
        ylims_caoi = (-5, 25)
    elif region == "NA":
        ylims_temp = (255, 310)
        ylims_caoi = (-15, 40)

    ylims_temp = (245, 310)
    ylims_caoi = (-10, 30)

    # PRESENT
    present_CAOI = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    present_SSTH = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    present_TH_850 = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)
    
    # FUTURE
    future_CAOI = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    future_SSTH = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    future_TH_850 = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)

    # Prepare data for boxplots
    data_temp = [
        present_SSTH.values.flatten(),
        future_SSTH.values.flatten(),
        present_TH_850.values.flatten(),
        future_TH_850.values.flatten(),
    ]
    data_caoi = [
        present_CAOI.values.flatten(),
        future_CAOI.values.flatten(),
    ]

    # filter out NaN values
    data_temp = [d[~np.isnan(d)] for d in data_temp]
    data_caoi = [d[~np.isnan(d)] for d in data_caoi]

    # Calculate mean differences (future - present)
    mean_diff_SSTH = np.nanmean(data_temp[1]) - np.nanmean(data_temp[0])
    mean_diff_TH850 = np.nanmean(data_temp[3]) - np.nanmean(data_temp[2])
    mean_diff_CAOI = np.nanmean(data_caoi[1]) - np.nanmean(data_caoi[0])

    # --- Plot both boxplots in the same figure ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 9), gridspec_kw={'width_ratios': [2, 1]})

    # SSTH & TH850
    bp_temp = ax1.boxplot(
        data_temp,
        positions=[0, 1, 2, 3],
        patch_artist=True,
        showfliers=False,
        widths=0.6,
    )

    ax1.set_ylabel(r"$\Theta$ (K)", fontsize=20)
    ax1.set_ylim(ylims_temp)
    ax1.tick_params(axis="both", labelsize=16)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)
    colors_temp = ["blue", "orange", "blue", "orange"]
    for patch, color in zip(bp_temp["boxes"], colors_temp):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax1.set_xticks([0.5, 2.5])
    xticklabels = [
        f"SST\nΔ={mean_diff_SSTH:.2f}K",
        f"850\nΔ={mean_diff_TH850:.2f}K"
    ]
    ax1.set_xticklabels(xticklabels, fontsize=20)
    print(ax1.get_ylim())

    # CAOI
    bp_caoi = ax2.boxplot(
        data_caoi,
        positions=[0, 1],
        patch_artist=True,
        showfliers=False,
        widths=0.6,
    )
    # Move ylabel to the right side
    ax2.set_ylabel("CAOI (K)", fontsize=20, labelpad=10)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax2.set_ylim(ylims_caoi)
    ax2.tick_params(axis="y", labelsize=16)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    colors_caoi = ["blue", "orange"]
    for patch, color in zip(bp_caoi["boxes"], colors_caoi):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax2.set_xticks([0.5])
    ax2.set_xticklabels([f"CAOI\nΔ={mean_diff_CAOI:.2f}K"], fontsize=20)
    print(ax2.get_ylim())

    plt.tight_layout()
    plt.savefig(f"{fig_path}/CAOI_contribution_combined_{region}.png", dpi=300, bbox_inches="tight", facecolor="white", transparent=False)
    plt.show()

def plot_caoi_contribution_djf(region):
    if region == "JS":
        ylims_temp = (245, 305)
        ylims_caoi = (-5, 25)
    elif region == "NA":
        ylims_temp = (255, 310)
        ylims_caoi = (-15, 40)

    ylims_temp = (245, 310)
    ylims_caoi = (-10, 26)

    # PRESENT
    present_CAOI = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    present_SSTH = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    present_TH_850 = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)
    
    # FUTURE
    future_CAOI = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    future_SSTH = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    future_TH_850 = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)

    # Prepare data for boxplots
    data_temp = [
        present_SSTH.values.flatten(),
        future_SSTH.values.flatten(),
        present_TH_850.values.flatten(),
        future_TH_850.values.flatten(),
    ]
    data_caoi = [
        present_CAOI.values.flatten(),
        future_CAOI.values.flatten(),
    ]

    # filter out NaN values
    data_temp = [d[~np.isnan(d)] for d in data_temp]
    data_caoi = [d[~np.isnan(d)] for d in data_caoi]

    # Calculate mean differences (future - present)
    mean_diff_SSTH = np.nanmean(data_temp[1]) - np.nanmean(data_temp[0])
    mean_diff_TH850 = np.nanmean(data_temp[3]) - np.nanmean(data_temp[2])
    mean_diff_CAOI = np.nanmean(data_caoi[1]) - np.nanmean(data_caoi[0])

    # --- Plot both boxplots in the same figure ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 9), gridspec_kw={'width_ratios': [2, 1]})

    # SSTH & TH850
    bp_temp = ax1.boxplot(
        data_temp,
        positions=[0, 1, 2, 3],
        patch_artist=True,
        showfliers=False,
        widths=0.6,
    )

    ax1.set_ylabel(r"$\Theta$ (K)", fontsize=20)
    # ax1.set_ylim(ylims_temp)
    ax1.tick_params(axis="both", labelsize=16)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)
    colors_temp = ["blue", "orange", "blue", "orange"]
    for patch, color in zip(bp_temp["boxes"], colors_temp):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax1.set_xticks([0.5, 2.5])
    xticklabels = [
        f"SST\nΔ={mean_diff_SSTH:.2f}K",
        f"850\nΔ={mean_diff_TH850:.2f}K"
    ]
    ax1.set_xticklabels(xticklabels, fontsize=20)
    print(ax1.get_ylim())

    # CAOI
    bp_caoi = ax2.boxplot(
        data_caoi,
        positions=[0, 1],
        patch_artist=True,
        showfliers=False,
        widths=0.6,
    )
    # Move ylabel to the right side
    ax2.set_ylabel("CAOI (K)", fontsize=20, labelpad=10)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax2.set_ylim(ylims_caoi)
    ax2.tick_params(axis="y", labelsize=16)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    colors_caoi = ["blue", "orange"]
    for patch, color in zip(bp_caoi["boxes"], colors_caoi):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax2.set_xticks([0.5])
    ax2.set_xticklabels([f"CAOI\nΔ={mean_diff_CAOI:.2f}K"], fontsize=20)
    print(ax2.get_ylim())

    plt.tight_layout()
    plt.savefig(f"{fig_path}/CAOI_contribution_combined_{region}_djf.png", dpi=300, bbox_inches="tight", facecolor="white", transparent=False)
    plt.show()

def plot_caoi_contribution_combined(region):
    if region == "JS":
        ylims_temp = (245, 305)
        ylims_caoi = (-5, 25)
    elif region == "NA":
        ylims_temp = (255, 310)
        ylims_caoi = (-15, 40)

    ylims_temp = (245, 310)
    ylims_caoi = (-10, 30)

    # EVENT DATA (foreground)
    event_CAOI_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    event_SSTH_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    event_TH_850_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)
    
    event_CAOI_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    event_SSTH_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    event_TH_850_future = get_event_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)

    # DJF DATA (background)
    djf_CAOI_present = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    djf_SSTH_present = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
    djf_TH_850_present = get_djf_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', member_ids)
    
    djf_CAOI_future = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'CAOI', member_ids)
    djf_SSTH_future = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', member_ids)
    djf_TH_850_future = get_djf_in_box_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', member_ids)

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

    # --- Plot both boxplots in the same figure ---
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
        medianprops={'linewidth': 1.5, 'alpha': 0.2},
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
        medianprops={'linewidth': 1.5},
    )
    colors_temp_event = ["blue", "orange", "blue", "orange"]
    for patch, color in zip(bp_temp_event["boxes"], colors_temp_event):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax1.set_ylabel(r"$\Theta$ (K)", fontsize=20)
    # ax1.set_ylim(ylims_temp)
    ax1.tick_params(axis="both", labelsize=16)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)
    ax1.set_xticks([0.5, 2.5])
    xticklabels = [
        f"SST\nΔ={mean_diff_SSTH:.2f}K",
        f"850\nΔ={mean_diff_TH850:.2f}K"
    ]
    ax1.set_xticklabels(xticklabels, fontsize=20)
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
        medianprops={'linewidth': 1.5, 'alpha': 0.2},
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
        medianprops={'linewidth': 1.5},
    )
    colors_caoi_event = ["blue", "orange"]
    for patch, color in zip(bp_caoi_event["boxes"], colors_caoi_event):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Move ylabel to the right side
    ax2.set_ylabel("CAOI (K)", fontsize=20, labelpad=10)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax2.set_ylim(ylims_caoi)
    ax2.tick_params(axis="y", labelsize=16)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    ax2.set_xticks([0.5])
    ax2.set_xticklabels([f"CAOI\nΔ={mean_diff_CAOI:.2f}K"], fontsize=20)
    print(ax2.get_ylim())

    # Add legend to distinguish between DJF and events
    legend_elements = [
        Patch(facecolor='gray', alpha=0.2, label='DJF'),
        Patch(facecolor='gray', alpha=0.7, label='Events')
    ]
    ax1.legend(handles=legend_elements, loc='upper left', fontsize=16)

    plt.tight_layout()
    plt.savefig(f"{fig_path}/CAOI_contribution_combined_overlay_{region}.png", dpi=300, bbox_inches="tight", facecolor="white", transparent=False)
    plt.show()


def plot_caoi_contribution_scatter(region):
    # data_present = get_percentile_in_box(data_dict[f'{region}_present'], region, 'present')
    # data_future = get_percentile_in_box(data_dict[f'{region}_future'], region, 'future')
    # present_SSTH, present_TH_850 = data_present.SSTH_masked, data_present.TH_850_masked
    # future_SSTH, future_TH_850 = data_future.SSTH_masked, data_future.TH_850_masked

    present_SSTH = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', event_mask_time, member_ids)
    present_TH_850 = get_masked_ensemble(data_dict[f'{region}_present'], region, 'present', 'TH_850_masked', event_mask_time, member_ids)
    future_SSTH = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'SSTH_masked', event_mask_time, member_ids)
    future_TH_850 = get_masked_ensemble(data_dict[f'{region}_future'], region, 'future', 'TH_850_masked', event_mask_time, member_ids)

    present_SSTH_filtered_mean = calculate_space_mean(present_SSTH.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3])))
    present_TH_850_filtered_mean = calculate_space_mean(present_TH_850.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3])))
    future_SSTH_filtered_mean = calculate_space_mean(future_SSTH.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3])))
    future_TH_850_filtered_mean = calculate_space_mean(future_TH_850.sel(lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3])))

    print(np.shape(present_SSTH_filtered_mean))
    print(np.shape(present_TH_850_filtered_mean))
    print(np.shape(future_SSTH_filtered_mean))
    print(np.shape(future_TH_850_filtered_mean))

    # Create scatter plot
    plt.figure(figsize=(8, 8), layout='constrained')
    plt.scatter(present_SSTH_filtered_mean, present_TH_850_filtered_mean, color='blue', label=f'present, n = {np.shape(present_SSTH_filtered_mean)[0]}', alpha=0.5)
    plt.plot(present_SSTH_filtered_mean.mean(skipna=True), present_TH_850_filtered_mean.mean(skipna=True), 'X', color='darkblue', markersize=10, label='mean')
    plt.scatter(future_SSTH_filtered_mean, future_TH_850_filtered_mean, color='orange', label=f'future, n = {np.shape(future_SSTH_filtered_mean)[0]}', alpha=0.5)
    plt.plot(future_SSTH_filtered_mean.mean(skipna=True), future_TH_850_filtered_mean.mean(skipna=True), 'X', color='darkorange', markersize=10, label='mean')

    plt.xlabel('$\Theta_{SST}$ (K)', fontsize=20)
    plt.ylabel('$\Theta_{850}$ (K)', fontsize=20)

    plt.xlim(280, 300)
    plt.xticks(fontsize=16)
    plt.ylim(271, 291)
    plt.yticks(fontsize=16)


    plt.legend(fontsize=20)
    plt.grid(alpha=0.5)
    plt.savefig(f"{fig_path}/CAOI_contribution_scatter_{region}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

def perform_eof_analysis(data_dict, variable, percentile_mask_dict, region, timeframe):
    """
    Perform EOF analysis on the data for a specific region and timeframe.

    Parameters:
        data_dict (dict): Dictionary containing the data for each region and timeframe.
        variable (str): Variable to analyze (e.g., "CAOI", "SST", "SSTH_masked", "TH_850_masked").
        percentile_mask_dict (dict): Dictionary containing the masks for 10th percentile for each region and timeframe.
        region (str): Region to analyze (e.g., "JS", "NA").
        timeframe (str): Timeframe to analyze (e.g., "present", "future").

    Returns:
        solver (Eof): EOF solver object containing the results of the analysis.
    """
    # Subset the data using the mask
    data = data_dict[f"{region}_{timeframe}"][variable]
    mask = percentile_mask_dict[f"{region}_{timeframe}"]
    data_subset = data.sel(time=mask)

    # Reshape the data for EOF analysis (time, space)
    # Flatten lat/lon into a single spatial dimension
    data_reshaped = data_subset.stack(space=("lat", "lon"))

    # Remove the mean along the time dimension (detrend)
    data_anomalies = data_reshaped - data_reshaped.mean(dim="time")

    # Perform EOF analysis
    solver = Eof(data_anomalies)

    return solver

def plot_eof_results(solver, data_dict, region, timeframe, variable, mode=1):
    """
    Plot the EOF pattern, explained variance, and principal components for a specific mode.

    Parameters:
        solver (Eof): EOF solver object.
        data_dict (dict): Dictionary containing the data for each region and timeframe.
        region (str): Region to analyze (e.g., "JS", "NA").
        timeframe (str): Timeframe to analyze (e.g., "present", "future").
        variable (str): Variable being analyzed (e.g., "CAOI", "SST").
        mode (int): EOF mode to plot (default is 1).
    """
    # Get the EOF pattern and explained variance
    eof_pattern = solver.eofsAsCovariance(neofs=mode).unstack("space")
    explained_variance = solver.varianceFraction(neigs=mode) * 100

    # Get the principal components (PCs)
    pcs = solver.pcs(npcs=mode)

    # Determine the color limits to center the colormap at 0
    vmin = eof_pattern.min().item()
    vmax = eof_pattern.max().item()
    abs_max = max(abs(vmin), abs(vmax))

    # Create a figure with two subplots: one for the EOF pattern and one for the PCs
    fig = plt.figure(figsize=(15, 18))

    # Create the first subplot with a specific projection
    ax1 = fig.add_subplot(2, 1, 1, projection=ccrs.PlateCarree())

    # Create the second subplot without a projection
    ax2 = fig.add_subplot(2, 1, 2)

    # Plot the EOF pattern
    ax = ax1
    ax.coastlines(color=(217/256, 217/256, 217/256))
    ax.set_facecolor("white")

    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 20, "color": "k"}
    gl.ylabel_style = {"size": 20, "color": "k"}
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(0, 100, 10))

    # Center the colormap at 0
    eof_pattern.plot(
        ax=ax,
        cmap="RdBu_r",
        robust=True,
        vmin=-abs_max,
        vmax=abs_max
    )
    ax.set_title(f"EOF Mode {mode} - {region} ({timeframe}) - {variable}\nExplained Variance: {explained_variance[0].item():.2f}%", fontsize=16)

    # Plot the PCs
    ax_pc = ax2
    pcs.plot(ax=ax_pc, linewidth=2)
    ax_pc.set_title(f"Principal Components (PCs) - Mode {mode}", fontsize=16)
    ax_pc.set_xlabel("Time", fontsize=14)
    ax_pc.set_ylabel("Amplitude", fontsize=14)
    ax_pc.grid(alpha=0.5)

    # Save the figure
    plt.savefig(f"{fig_path}/EOF_and_PCs_{region}_{timeframe}_{variable}_mode{mode}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

def plot_histogram_event(region):
    durations_2000_list = []
    durations_2090_list = []
    n_2000 = 0
    n_2090 = 0

    # Loop over all ensemble members
    for member_id in member_ids:
        file_events = f'/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/events_d/dates_CAO_max_intens_1990:2100_{region}'
        events = pd.read_csv(file_events, sep=' ', header=None)
        events.columns = ['time', 'CAOI_ano', 'event_ID', 'duration']
        events['time'] = pd.to_datetime(events['time'], format='%Y-%m-%d-%H')
        events['year'] = events['time'].dt.year

        durations_2000 = events[events["year"] <= 2000]["duration"]
        durations_2090 = events[events["year"] >= 2090]["duration"]

        n_2000 += durations_2000.shape[0]
        n_2090 += durations_2090.shape[0]
        durations_2000_list.append(durations_2000)
        durations_2090_list.append(durations_2090)

        #print duration and 'time' of longest event for each member
        if not durations_2000.empty:
            max_duration_2000 = durations_2000.max()
            max_event_2000 = events[(events["year"] <= 2000) & (events["duration"] == max_duration_2000)]
            print(f"Member {member_id} - Longest event (<=2000): Duration = {max_duration_2000} hours, Time = {max_event_2000['time'].values[0]}")
        if not durations_2090.empty:
            max_duration_2090 = durations_2090.max()
            max_event_2090 = events[(events["year"] >= 2090) & (events["duration"] == max_duration_2090)]
            print(f"Member {member_id} - Longest event (>=2090): Duration = {max_duration_2090} hours, Time = {max_event_2090['time'].values[0]}")

    # Stack all durations from all members
    all_durations_2000 = np.concatenate([d.values for d in durations_2000_list])
    all_durations_2090 = np.concatenate([d.values for d in durations_2090_list])

    mean_duration_2000 = np.mean(all_durations_2000)
    mean_duration_2090 = np.mean(all_durations_2090)
    print(f"Mean event duration (<=2000): {mean_duration_2000:.2f} hours")
    print(f"Mean event duration (>=2090): {mean_duration_2090:.2f} hours")
    
    min_duration = min(all_durations_2000.min(), all_durations_2090.min())
    max_duration = max(all_durations_2000.max(), all_durations_2090.max())


    # Create bins with width 24 hours
    # Start from the first multiple of 24 at or below min_duration
    bin_start = int((min_duration // 24) * 24)
    # End at the first multiple of 24 above max_duration
    bin_end = int(((max_duration // 24) + 1) * 24)
    bins = np.arange(bin_start, bin_end + 24, 24)
    
    # Create the plot with more spacious layout
    plt.figure(figsize=(10, 6), layout='constrained')
    
    n_hist_2000, _, _ = plt.hist(all_durations_2000, bins=bins, alpha=0)  # Get counts without plotting
    n_hist_2090, _, _ = plt.hist(all_durations_2090, bins=bins, alpha=0)  # Get counts without plotting
    plt.clf()  # Clear the invisible histograms
    
    # Plot histograms side by side with more spacing
    width = (bins[1] - bins[0]) * 0.35  # Make bars slightly narrower for more spacing
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    # Normalize the counts to get relative frequencies
    n_hist_2000_norm = n_hist_2000 / n_2000 if n_2000 > 0 else n_hist_2000
    n_hist_2090_norm = n_hist_2090 / n_2090 if n_2090 > 0 else n_hist_2090
    
    bars_2000 = plt.bar(bin_centers - width/2, n_hist_2000_norm, width=width, 
                       label=f'present (n = {n_2000})', color='blue', alpha=0.7)
    bars_2090 = plt.bar(bin_centers + width/2, n_hist_2090_norm, width=width, 
                       label=f'future (n = {n_2090})', color='orange', alpha=0.7)

    # Add count labels above bars
    for i, (bar_2000, bar_2090, count_2000, count_2090) in enumerate(zip(bars_2000, bars_2090, n_hist_2000, n_hist_2090)):
        # Add count for 2000 data if non-zero
        if count_2000 > 0:
            x = bar_2000.get_x() + bar_2000.get_width() / 2
            plt.text(x, bar_2000.get_height(), str(int(count_2000)), ha='center', va='bottom', fontsize=10)
        
        # Add count for 2090 data if non-zero
        if count_2090 > 0:
            x = bar_2090.get_x() + bar_2090.get_width() / 2
            plt.text(x, bar_2090.get_height(), str(int(count_2090)), ha='center', va='bottom', fontsize=10)
    
    # Customize the plot with matching style
    plt.xlabel("duration [h]", fontsize=20)
    plt.ylabel("frequency", fontsize=20)
    plt.legend(fontsize=20)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    
    # Set x-axis ticks to show the bin boundaries clearly
    plt.xticks(bins, [f'{int(x)}' for x in bins], fontsize=16, rotation=45)
    plt.yticks(fontsize=16)
    
    plt.savefig(f"{fig_path}/event_duration_histogram_{region}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

def plot_CAOI_difference(region):
    CAOI_present = get_masked_ensemble(data_dict[f'{region}_present'], region, "present", 'CAOI', event_mask_time, member_ids)
    CAOI_present_corrected = get_masked_ensemble(data_dict[f'{region}_present'], region, "present", 'CAOI', event_mask_time_corrected, member_ids)
    CAOI_diff = calculate_time_mean(CAOI_present_corrected) - calculate_time_mean(CAOI_present)

    # Use the same plotting setup as plot_single_map to visualize the difference
    diff_data = CAOI_diff
    lons = CAOI_present.lon
    lats = CAOI_present.lat
    region_extent = map_extent[region]

    fig, ax = plt.subplots(
        1, 1,
        figsize=get_figsize_from_map_extent(region_extent, base_height=7),
        subplot_kw=dict(projection=ccrs.PlateCarree()),
        layout='constrained'
    )
    ax.set_extent(region_extent, crs=ccrs.PlateCarree())
    ax.coastlines(color='black')
    ax.set_facecolor("white")
    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 20, "color": "k"}
    gl.ylabel_style = {"size": 20, "color": "k"}
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(0, 100, 10))

    cmap, n_main = cmap_pcolormesh('RdBu_r', np.arange(-3, 3.5, 0.5))
    norm = BoundaryNorm(np.arange(-3, 3.5, 0.5), ncolors=n_main, clip=False)
    im = ax.pcolormesh(lons, lats, diff_data, shading='auto',
                       transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)

    divider = make_axes_locatable(ax)
    shrink = get_colorbar_shrink(region_extent)
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', shrink=shrink, pad=0.05, extend='both')
    cbar.set_label('CAOI difference (corrected - original) (K)', size=20)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_ticks(np.arange(-3, 4, 1))

    plt.title('Difference between corrected and original CAOI masking', fontsize=18)
    plt.savefig(f"{fig_path}/CAOI_masking_difference_{region}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()
    plt.close()

# #%%
# write_data("JS", "present")
# write_data("JS", "future")
# write_data("NA", "present")
# write_data("NA", "future")
#%%
## read data (all ensemble members)
for region in regions:
    for timeframe in timeframes:
        ds_ens = read_ensemble_data(region, timeframe, member_ids)

#%%
quantiles_dict = {
    "NA_present": calculate_quantiles('NA', 'present'),
    "NA_future": calculate_quantiles('NA', 'future'),
    "JS_present": calculate_quantiles('JS', 'present'),
    "JS_future": calculate_quantiles('JS', 'future'),
    # add other regions/timeframes as needed
}
#%%
data_dict = {
    "NA_present": read_ensemble_data("NA", "present", member_ids),
    "NA_future": read_ensemble_data("NA", "future", member_ids),
    "JS_present": read_ensemble_data("JS", "present", member_ids),
    "JS_future": read_ensemble_data("JS", "future", member_ids)
    # add other regions/timeframes as needed
}

#%%

data_dict_precip = {
    "NA_present": read_ensemble_data_precip("NA", "present", member_ids),
    "NA_future": read_ensemble_data_precip("NA", "future", member_ids),
    "JS_present": read_ensemble_data_precip("JS", "present", member_ids),
    "JS_future": read_ensemble_data_precip("JS", "future", member_ids)
    # add other regions/timeframes as needed
}

#%%

data_dict_ice = {
    "NA_present": read_ensemble_data_ice("NA", "present", member_ids),
    "NA_future": read_ensemble_data_ice("NA", "future", member_ids),
    "JS_present": read_ensemble_data_ice("JS", "present", member_ids),
    "JS_future": read_ensemble_data_ice("JS", "future", member_ids)
    # add other regions/timeframes as needed
}

#%% 
# TODO's


#%%
# Define regions and timeframes
regions = ["JS", "NA"]
timeframes = ["present", "future"]
variables = ["CAOI", "SSTH_masked", "TH_850_masked"]

# # EOF analysis
# for region in regions:
#     for variable in variables:
#         for timeframe in timeframes:
#             # Perform EOF analysis
#             print(f"Performing EOF analysis for {region} during {timeframe} for variable {variable}")

#             # Perform EOF analysis
#             solver = perform_eof_analysis(data_dict, variable, percentile_mask_dict, region, timeframe)

#             # Plot EOF results for the first mode
#             plot_eof_results(solver, data_dict, region, timeframe, variable, mode=1)


#%%
## good plots
plot_map_10th("JS")
plot_map_10th("NA")
plot_map_frequency("JS")
plot_map_frequency("NA")
plot_map_frequency_10th("JS")
plot_map_frequency_10th("NA")
plot_map_fluxes_10th("JS")
plot_map_fluxes_10th("NA")
plot_map_djf("JS")
plot_map_djf("NA")
plot_boxplot("JS")
plot_boxplot("NA")
plot_histogram("JS")
plot_histogram("NA")
plot_caoi_contribution_bar("JS")
plot_caoi_contribution_bar("NA")
plot_caoi_contribution("JS")
plot_caoi_contribution("NA")
plot_caoi_contribution_scatter("JS")
plot_caoi_contribution_scatter("NA")
plot_histogram_event("JS")
plot_histogram_event("NA")


#%%
regions = ["JS", "NA"]
varnames = ['EP', 'PRECIP', 'EVAP']
for region in regions:
    for varname in varnames:
        for lag in [0, 1, 2]:
            print(f"Plotting {varname} for {region} with lag {lag}")
            plot_map_10th(region, varname=varname, lag=lag)

#%%
regions = ["JS", "NA"]
for region in regions:
    plot_map_10th(region, varname='CAOI', lag=0)

# %%

# Define extents
extent_JS = [110, 160, 20, 85]
extent_NA = [260, 310, 20, 85]

central_longitudes = {
    "JS": 140,  # for map projection
    "NA": 290   # for map projection
}
def create_map_projection(region, central_longitudes):
    """Create map projection for given region"""
    map_projection = ccrs.LambertConformal(central_longitude=central_longitudes[region], central_latitude=75)
    return map_projection

# Use the same facecolor as in plot_single_map
land_facecolor = "white"

def plot_region_box(extent, box_extent, color, label):
    fig = plt.figure(figsize=get_figsize_from_map_extent(extent, base_height=10), layout='constrained')
    ax = plt.axes(projection=create_map_projection(label, central_longitudes))
    ax.set_facecolor("white")

    gl = ax.gridlines(draw_labels=True, dms=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 20, "color": "k"}
    gl.ylabel_style = {"size": 20, "color": "k"}
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(0, 80, 10))
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor=(217/256, 217/256, 217/256))

    if label == 'NA':
        # Adjust box for NA region crossing the 0° meridian
        box_extent = [box_extent[0]-360, box_extent[1]-360, box_extent[2], box_extent[3]]

    # Draw the region box
    rect = Rectangle(
        (box_extent[0], box_extent[2]),
        box_extent[1] - box_extent[0],
        box_extent[3] - box_extent[2],
        linewidth=2, edgecolor=color, facecolor='none', label=label, transform=ccrs.PlateCarree()
    )
    ax.add_patch(rect)

    lons = data_dict_ice[f"{label}_present"].sel(lon=slice(box[label][0], box[label][1]), lat=slice(box[label][2], box[label][3])).lon.values
    lats = data_dict_ice[f"{label}_present"].sel(lon=slice(box[label][0], box[label][1]), lat=slice(box[label][2], box[label][3])).lat.values
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    ax.scatter(lon_grid, lat_grid, color=color, s=1, transform=ccrs.PlateCarree(), label=f"{label} grid points")
    plt.savefig(f"{fig_path}/region_box_{label}.png", dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
    plt.show()

# JS map
plot_region_box(extent_JS, box["JS"], 'darkblue', 'JS')

# NA map
plot_region_box(extent_NA, box["NA"], 'darkblue', 'NA')


#%%
def check_sea_ice(varname, region, timeframe_name, lag):
    if varname in ['PRECIP', 'EVAP', 'EP']:
        # For precipitation, evaporation, and P-E, use the map extent for the box
        if region == "NA":
            # For North Atlantic, use a different projection
            map_projection = ccrs.AzimuthalEquidistant(central_longitude=-40, central_latitude=50)
            extent = [map_extent_precip[region][0]-360, map_extent_precip[region][1]-15, map_extent_precip[region][2], map_extent_precip[region][3]]
        elif region == 'JS':
            map_projection = ccrs.AzimuthalEquidistant(central_longitude=190, central_latitude=50)
            extent = map_extent_precip[region]
        
        if lag == 0:
            event_mask_time_func = event_plus0_mask_time
        elif lag == 1:
            event_mask_time_func = event_plus1_mask_time
        elif lag == 2:
            event_mask_time_func = event_plus2_mask_time
        area_slice = read_extent_precip[region] # to match sea ice data with precipitation data
    else:
        # For other variables, use the box defined in the VAR_CONFIG
        map_projection = ccrs.PlateCarree()
        extent = map_extent[region]
        event_mask_time_func = event_mask_time
        area_slice = map_extent[region] # to match sea ice data with CAOI data

    if timeframe_name == "10th_present":
        ds_name = 'present'
        # SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP']:
            SLP = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SLP_ano", event_mask_time_func, member_ids)
            time_slice = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
        else:
            time_slice = djf_mask(data_dict_ice[f"{region}_present"], 'present')
        
        print(f"Using time slice for {region} {ds_name}: {time_slice}")
        
        # Apply the same longitude selection logic as for precipitation data
        sea_ice_full = data_dict_ice[f"{region}_present"].sel(
            lat=slice(area_slice[2], area_slice[3]),
            time=time_slice
        )
        
        # Handle longitude selection for regions that cross 0° meridian
        if area_slice[1] < area_slice[0]:  # Longitude wraps around (like NA region)
            # Select west part (e.g., 275° to 360°) and east part (e.g., 0° to 15°)
            sea_ice_west = sea_ice_full.sel(lon=slice(area_slice[0], 360))
            sea_ice_east = sea_ice_full.sel(lon=slice(0, area_slice[1]))
            sea_ice_sliced = xr.concat([sea_ice_west, sea_ice_east], dim='lon')
        else:
            # Normal case - no longitude wrapping
            sea_ice_sliced = sea_ice_full.sel(lon=slice(area_slice[0], area_slice[1]))
        
        sea_ice = get_masked_ensemble(sea_ice_sliced, region, 'present', "ICEFRAC", event_mask_time_func, member_ids)
        print(f"Using sea ice data for {region} {ds_name}: {sea_ice}")
        print(calculate_time_mean(sea_ice).mean())

    elif timeframe_name == "10th_future":
        ds_name = 'future'
        # SST = get_masked_ensemble(data_dict[f"{region}_{ds_name}"], region, ds_name, "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP']:
            SLP = get_masked_ensemble(data_dict_precip[f"{region}_{ds_name}"], region, ds_name, "SLP_ano", event_mask_time_func, member_ids)
            time_slice = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
        else:
            time_slice = djf_mask(data_dict_ice[f"{region}_{ds_name}"], ds_name)
        
        print(f"Using time slice for {region} {ds_name}: {time_slice}")

        # Apply the same longitude selection logic as for precipitation data
        sea_ice_full = data_dict_ice[f"{region}_{ds_name}"].sel(
            lat=slice(area_slice[2], area_slice[3]),
            time=time_slice
        )
        
        # Handle longitude selection for regions that cross 0° meridian
        if area_slice[1] < area_slice[0]:  # Longitude wraps around (like NA region)
            # Select west part (e.g., 275° to 360°) and east part (e.g., 0° to 15°)
            sea_ice_west = sea_ice_full.sel(lon=slice(area_slice[0], 360))
            sea_ice_east = sea_ice_full.sel(lon=slice(0, area_slice[1]))
            sea_ice_sliced = xr.concat([sea_ice_west, sea_ice_east], dim='lon')
        else:
            # Normal case - no longitude wrapping
            sea_ice_sliced = sea_ice_full.sel(lon=slice(area_slice[0], area_slice[1]))
        
        sea_ice = get_masked_ensemble(sea_ice_sliced, region, ds_name, "ICEFRAC", event_mask_time_func, member_ids)
        print(f"Using sea ice data for {region} {ds_name}: {sea_ice}")
        print(calculate_time_mean(sea_ice).mean())

    elif timeframe_name == "10th_diff":
        # SST_present = get_masked_ensemble(data_dict[f"{region}_present"], region, "present", "SST", event_mask_time_func, member_ids)
        # SST_future = get_masked_ensemble(data_dict[f"{region}_future"], region, "future", "SST", event_mask_time_func, member_ids)
        if varname in ['PRECIP', 'EVAP', 'EP']:
            SLP_present = get_masked_ensemble(data_dict_precip[f"{region}_present"], region, "present", "SLP_ano", event_mask_time_func, member_ids)
            SLP_future = get_masked_ensemble(data_dict_precip[f"{region}_future"], region, "future", "SLP_ano", event_mask_time_func, member_ids)
            time_slice_present = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
            time_slice_future = slice(None, None)  # Use all times for precipitation, evaporation, and P-E
        else:
            time_slice_present = djf_mask(data_dict_ice[f"{region}_present"], "present")
            time_slice_future = djf_mask(data_dict_ice[f"{region}_future"], "future")

        # Handle present data
        sea_ice_full_present = data_dict_ice[f"{region}_present"].sel(
            lat=slice(area_slice[2], area_slice[3]),
            time=time_slice_present
        )
        
        if area_slice[1] < area_slice[0]:  # Longitude wraps around
            sea_ice_west_present = sea_ice_full_present.sel(lon=slice(area_slice[0], 360))
            sea_ice_east_present = sea_ice_full_present.sel(lon=slice(0, area_slice[1]))
            sea_ice_sliced_present = xr.concat([sea_ice_west_present, sea_ice_east_present], dim='lon')
        else:
            sea_ice_sliced_present = sea_ice_full_present.sel(lon=slice(area_slice[0], area_slice[1]))
        
        # Handle future data
        sea_ice_full_future = data_dict_ice[f"{region}_future"].sel(
            lat=slice(area_slice[2], area_slice[3]),
            time=time_slice_future
        )
        
        if area_slice[1] < area_slice[0]:  # Longitude wraps around
            sea_ice_west_future = sea_ice_full_future.sel(lon=slice(area_slice[0], 360))
            sea_ice_east_future = sea_ice_full_future.sel(lon=slice(0, area_slice[1]))
            sea_ice_sliced_future = xr.concat([sea_ice_west_future, sea_ice_east_future], dim='lon')
        else:
            sea_ice_sliced_future = sea_ice_full_future.sel(lon=slice(area_slice[0], area_slice[1]))
        
        sea_ice_present = get_masked_ensemble(sea_ice_sliced_present, region, "present", "ICEFRAC", event_mask_time_func, member_ids)
        sea_ice_future = get_masked_ensemble(sea_ice_sliced_future, region, "future", "ICEFRAC", event_mask_time_func, member_ids)
    calculate_time_mean(sea_ice).plot()

#%% test boxplot
def plot_caoi_contribution_combined(region):
    if region == "JS":
        ylims_temp = (245, 305)
        ylims_caoi = (-5, 25)
    elif region == "NA":
        ylims_temp = (255, 310)
        ylims_caoi = (-15, 40)

    ylims_temp = (245, 310)
    ylims_caoi = (-10, 30)

    # EVENT DATA (foreground)
    event_CAOI_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'CAOI', member_ids)
    event_CAOI_present = calculate_space_mean(event_CAOI_present)
    event_SSTH_present = get_event_in_box_ensemble(data_dict[f'{region}_present'], region, 'present', 'SSTH_masked', member_ids)
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

    # --- Plot both boxplots in the same figure ---
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
        medianprops={'linewidth': 1.5, 'alpha': 0.2},
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
        medianprops={'linewidth': 1.5},
    )
    colors_temp_event = ["blue", "orange", "blue", "orange"]
    for patch, color in zip(bp_temp_event["boxes"], colors_temp_event):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax1.set_ylabel(r"$\Theta$ [K]", fontsize=20)
    # ax1.set_ylim(ylims_temp)
    ax1.tick_params(axis="both", labelsize=16)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)
    ax1.set_xticks([0.5, 2.5])
    xticklabels = [
        f"SST\nΔ={mean_diff_SSTH:.2f}K",
        f"850\nΔ={mean_diff_TH850:.2f}K"
    ]
    ax1.set_xticklabels(xticklabels, fontsize=20)
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
        medianprops={'linewidth': 1.5, 'alpha': 0.2},
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
        medianprops={'linewidth': 1.5},
    )
    colors_caoi_event = ["blue", "orange"]
    for patch, color in zip(bp_caoi_event["boxes"], colors_caoi_event):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Move ylabel to the right side
    ax2.set_ylabel("CAOI [K]", fontsize=20, labelpad=10)
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax2.set_ylim(ylims_caoi)
    ax2.tick_params(axis="y", labelsize=16)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    ax2.set_xticks([0.5])
    ax2.set_xticklabels([f"CAOI\nΔ={mean_diff_CAOI:.2f}K"], fontsize=20)
    print(ax2.get_ylim())

    # Add legend to distinguish between DJF and events
    legend_elements = [
        Patch(facecolor='gray', alpha=0.2, label='DJF'),
        Patch(facecolor='gray', alpha=0.7, label='Events'),
        Patch(facecolor='blue', alpha=0.7, label='present'),
        Patch(facecolor='orange', alpha=0.7, label='future')
    ]
    ax1.legend(handles=legend_elements, loc='upper left', fontsize=16)

    plt.tight_layout()
    plt.savefig(f"{fig_path}/CAOI_contribution_combined_overlay_{region}.png", dpi=300, bbox_inches="tight", facecolor="white", transparent=False)
    plt.show()
# %% box vs map

##CAOI
#map
tcao_map_present = get_masked_ensemble(data_dict["JS_present"], "JS", "present", "CAOI", event_mask_time, member_ids)
tcao_map_mean_present = calculate_time_mean(tcao_map_present)
tcao_box_mean_present = tcao_map_mean_present.sel(lon=slice(box["JS"][0], box["JS"][1]), lat=slice(box["JS"][2], box["JS"][3]))
tcao_mean_box_mean_present = calculate_space_mean(tcao_box_mean_present)
#9.62037744

#boxplot
tcao_box_present = get_event_in_box_ensemble(data_dict["JS_present"], "JS", "present", "CAOI", member_ids)
tcao_box_mean_present = calculate_space_mean(tcao_box_present)
tcao_mean_box_mean_present = tcao_box_mean_present.mean()
#9.88975352

#SSTH_masked
#map
tcao_map_present = get_masked_ensemble(data_dict["JS_present"], "JS", "present", "SSTH_masked", event_mask_time, member_ids)
tcao_map_mean_present = calculate_time_mean(tcao_map_present)
tcao_box_mean_present = tcao_map_mean_present.sel(lon=slice(box["JS"][0], box["JS"][1]), lat=slice(box["JS"][2], box["JS"][3]))
tcao_mean_box_mean_present = calculate_space_mean(tcao_box_mean_present)
#283.21883564

#boxplot
tcao_box_present = get_event_in_box_ensemble(data_dict["JS_present"], "JS", "present", "SSTH_masked", member_ids)
tcao_box_mean_present = calculate_space_mean(tcao_box_present)
tcao_mean_box_mean_present = tcao_box_mean_present.mean()
#284.99171036

##TH_850_masked
#map
tcao_map_present = get_masked_ensemble(data_dict["JS_present"], "JS", "present", "TH_850_masked", event_mask_time, member_ids)
tcao_map_mean_present = calculate_time_mean(tcao_map_present)
tcao_box_mean_present = tcao_map_mean_present.sel(lon=slice(box["JS"][0], box["JS"][1]), lat=slice(box["JS"][2], box["JS"][3]))
tcao_mean_box_mean_present = calculate_space_mean(tcao_box_mean_present)
#273.59845821

#boxplot
tcao_box_present = get_event_in_box_ensemble(data_dict["JS_present"], "JS", "present", "TH_850_masked", member_ids)
tcao_box_mean_present = calculate_space_mean(tcao_box_present)
tcao_mean_box_mean_present = tcao_box_mean_present.mean()
#275.10195684

##CAOI
#map
tcao_map_future = get_masked_ensemble(data_dict["JS_future"], "JS", "future", "CAOI", event_mask_time, member_ids)
tcao_map_mean_future = calculate_time_mean(tcao_map_future)
tcao_box_mean_future = tcao_map_mean_future.sel(lon=slice(box["JS"][0], box["JS"][1]), lat=slice(box["JS"][2], box["JS"][3]))
tcao_mean_box_mean_future = calculate_space_mean(tcao_box_mean_future)
#8.28213165

#boxplot
tcao_box_future = get_event_in_box_ensemble(data_dict["JS_future"], "JS", "future", "CAOI", member_ids)
tcao_box_mean_future = calculate_space_mean(tcao_box_future)
tcao_mean_box_mean_future = tcao_box_mean_future.mean()
#8.2748464

#%% signficance
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

def debug_significance_data(present_data, future_data):
    """Debug function to check data before significance testing"""
    print("=== Data Debugging ===")
    print(f"Present data shape: {present_data.shape}")
    print(f"Future data shape: {future_data.shape}")
    print(f"Present data type: {present_data.dtype}")
    print(f"Future data type: {future_data.dtype}")
    
    # Check for NaN values
    present_nans = np.isnan(present_data).sum()
    future_nans = np.isnan(future_data).sum()
    print(f"Present NaNs: {present_nans.item()} / {present_data.size}")
    print(f"Future NaNs: {future_nans.item()} / {future_data.size}")
    
    # Check data ranges
    present_valid = present_data[~np.isnan(present_data)]
    future_valid = future_data[~np.isnan(future_data)]
    
    if len(present_valid) > 0:
        print(f"Present data range: {present_valid.min().item():.3f} to {present_valid.max().item():.3f}")
        print(f"Present data mean: {present_valid.mean().item():.3f}")
    
    if len(future_valid) > 0:
        print(f"Future data range: {future_valid.min().item():.3f} to {future_valid.max().item():.3f}")
        print(f"Future data mean: {future_valid.mean().item():.3f}")

# %%
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

# %%
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

    # Plot present (no stippling)
    plot_single_map(
        data_2d=CAOI_present_freq,
        lons=lons, lats=lats, region=region,
        timeframe_name='present', extent=box,
        cmap_name='bilbao', boundaries=np.arange(0, 110, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='frequency exceedance 8K [%]',
        extend='neither', save_path=True,
        significance_mask=None  # No stippling for present
    )

    # Plot future (no stippling)
    plot_single_map(
        data_2d=CAOI_future_freq,
        lons=lons, lats=lats, region=region,
        timeframe_name='future', extent=box,
        cmap_name='bilbao', boundaries=np.arange(0, 110, 10),
        varname='CAOI_8Kfrequency_10th',
        labname='frequency exceedance 8K [%]',
        extend='neither', save_path=True,
        significance_mask=None  # No stippling for future
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
            f"$\Theta_{{{SST_text}}}$\nΔ={mean_diff_SSTH:.2f}K{asterisk_SSTH}",
            f"$\Theta_{{{TH_text}}}$\nΔ={mean_diff_TH850:.2f}K{asterisk_TH850}"
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
# %%

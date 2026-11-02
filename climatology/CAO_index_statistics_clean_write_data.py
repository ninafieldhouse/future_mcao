# %%
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from lagranto.plotting import plot_trajs
from lagranto import Tra 
import matplotlib.pyplot as plt
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
from matplotlib.lines import Line2D  # Make sure this import is included
import xarray as xr
import matplotlib.dates as mdates  # Import mdates for date formatting
import cftime
import os
from cartopy.io import shapereader
from cartopy.feature import LAND
import cftime  # Add this import
from eofs.xarray import Eof
from matplotlib.patches import Rectangle
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
import config

# %%
member_ids = [1000, 1100, 1200, 1300, 1400]
regions = ["JS", "NA"]
timeframes = ["present", "future"]
base_path = f"{config.DATA_ROOT}/cao/clim/"
fig_path = f"{config.PLOTS_ROOT}/statistics"

map_extent = {
    "JS": [120, 160, 20, 70],  # for map
    "NA": [275, 305, 20, 55]   # for map
}

map_extent = {
    "JS": [120, 160, 20, 65],  # for map
    "NA": [270, 310, 20, 65]   # for map
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
def read_masked_CAOI_data(region, timeframe):
    """Read CAOI data masked by sea-ice, land, box (time, lat, lon)"""
    path = f"{base_path}/{region}/"
    ds_C_ano = xr.open_dataset(path + f"yearly_hour_anomalies_{timeframe}", decode_times=True)
    ds_C = xr.open_dataset(path + f"all_times_{timeframe}", decode_times=True)

    CAOI = -ds_C.CAOI[:, 0, :, :]
    CAOI_djf = -ds_C.sel(time=djf_mask(ds_C, timeframe)).CAOI[:, 0, :, :]
    # CAOI anomalies masked by sea-ice, land, box (time, lat, lon)
    CAOI_ano = ds_C_ano.CAOI[:, 0, :, :]
    CAOI_djf_ano = ds_C_ano.sel(time=djf_mask(ds_C_ano, timeframe)).CAOI[:, 0, :, :]

    return CAOI_djf, CAOI_djf_ano

def write_data(region, timeframe, member_id):
    """Write DJF netcdf data masked by sea-ice, land, NOT BOX (time, lat, lon)"""
    if timeframe == "present":
        years = range(1990, 2001)  # for example: 1990 to 2000
        cesm_type = 'HIST'
    elif timeframe == "future":
        years = range(2090, 2101)
        cesm_type = 'SSP370'

    cesm = f"{config.CESM_ARCHIVE}/b.e212.B{cesm_type}cmip6.f09_g17.{member_id}/archive/atm/hist/"

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

def write_data_time_extended(region, timeframe, member_id):
    """Write netcdf data masked by sea-ice, land, NOT BOX (time, lat, lon)"""
    if timeframe == "present":
        years = range(1990, 2001)  # for example: 1990 to 2000
        cesm_type = 'HIST'
    elif timeframe == "future":
        years = range(2090, 2101)
        cesm_type = 'SSP370'

    cesm = f"{config.CESM_ARCHIVE}/b.e212.B{cesm_type}cmip6.f09_g17.{member_id}/archive/atm/hist/"

    print('Opening files')
    # Construct file paths for each year
    files_to_open = [
        os.path.join(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.{member_id}.cam.h3.{year}-01-01-00000.nc")
        for year in years
    ]
    print('Opened files')

    # Open and combine the datasets
    ds = xr.open_mfdataset(files_to_open, combine='by_coords', decode_times=True)
    SST = ds.sel(time=ds['time'].dt.month.isin([11, 12, 1, 2, 3]), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).SST
    print('Opened SST')
    T = ds.sel(time=ds['time'].dt.month.isin([11, 12, 1, 2, 3]), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).T[:,25,:,:]
    print('Opened T')
    PS = ds.sel(time=ds['time'].dt.month.isin([11, 12, 1, 2, 3]), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).PS
    print('Opened PS')
    SHFLX = ds.sel(time=ds['time'].dt.month.isin([11, 12, 1, 2, 3]), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).SHFLX
    print('Opened SHFLX')
    LHFLX = ds.sel(time=ds['time'].dt.month.isin([11, 12, 1, 2, 3]), lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])).LHFLX
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
    file_path = f"{output_path}data_time_extended_{timeframe}"

    # Check if file exists and if you have write permissions
    if os.path.exists(file_path):
        if os.access(file_path, os.W_OK):
            # Remove the file before overwriting (to avoid permission errors)
            os.remove(file_path)
        else:
            raise PermissionError(f"No write permission for {file_path}")
    ds_to_save.to_netcdf(file_path, mode='w', format='NETCDF4', engine='netcdf4')

def write_data_ice(region, timeframe, member_id):
    if timeframe == "present":
        years = range(1990, 2001)  # for example: 1990 to 2000
        cesm_type = 'HIST'
    elif timeframe == "future":
        years = range(2090, 2101)
        cesm_type = 'SSP370'

    cesm = f"{config.CESM_ARCHIVE}/b.e212.B{cesm_type}cmip6.f09_g17.{member_id}/archive/atm/hist/"

    print('Opening files')
    # Construct file paths for each year
    files_to_open = [
        os.path.join(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.{member_id}.cam.h4.{year}-01-01-00000.nc")
        for year in years
    ]
    print('Opened files')

    # Open and combine the datasets
    ds = xr.open_mfdataset(files_to_open, combine='by_coords', decode_times=True)
    ICEFRAC = ds.sel(time=ds['time'].dt.month.isin([11, 12, 1, 2, 3])
                #  , lon=slice(map_extent[region][0], map_extent[region][1]), lat=slice(map_extent[region][2], map_extent[region][3])
                 ).ICEFRAC
    print('Opened ICEFRAC')

    ds_to_save = xr.Dataset(
            {
                "ICEFRAC": (["time", "lat", "lon"], ICEFRAC.data),
            },
            coords={
                "time": ICEFRAC.time,
                "lat": ICEFRAC.lat,
                "lon": ICEFRAC.lon,
            },
        )
    print('Saving data')

    # Save the dataset to a netCDF file
    output_path = f"{base_path}/{region}/"
    os.makedirs(output_path, exist_ok=True)
    file_path = f"{output_path}data_ice_{timeframe}"

    # Check if file exists and if you have write permissions
    if os.path.exists(file_path):
        if os.access(file_path, os.W_OK):
            # Remove the file before overwriting (to avoid permission errors)
            os.remove(file_path)
        else:
            raise PermissionError(f"No write permission for {file_path}")
    ds_to_save.to_netcdf(file_path, mode='w', format='NETCDF4', engine='netcdf4')


for member_id in member_ids:
    base_path = f"{config.DATA_ROOT}/{member_id}/cao/clim"
    for region in regions:
        for timeframe in timeframes:
            print(f"Processing {region} {timeframe} for member {member_id}")
            write_data(region, timeframe, member_id)
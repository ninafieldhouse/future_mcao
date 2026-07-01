# %%
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import netCDF4 as nc
from dypy.small_tools import interpolate
from geopy.distance import great_circle
import metpy.calc as mpcalc
from metpy.units import units
import xarray as xr
import itertools
import sys
from numba import jit
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import datetime as dt

@jit(nopython=True, boundscheck=False, cache=True, fastmath=True)
def v_interpolate_loop(z_out, data_out, z_in, data_in):
    ny, nx = z_out.shape
    for i in range(nx):
        for j in range(ny):
            data_out[j, i] = np.interp(
                z_out[j, i],
                z_in[:, j, i],
                data_in[:, j, i],
            )

def vertical_interpolate_3d(z_out, z_in, data_in):
    ny, nx = z_out.shape
    data_out = np.zeros((ny, nx))
    v_interpolate_loop(z_out, data_out, z_in, data_in)
    return data_out

# %%
# -----------------------------------
# define boundary conditions
# -----------------------------------
member_id = sys.argv[5]
# member_id = 1400
region = sys.argv[3]
# region = 'NA'
if region == "NA":
    box = [285,295,30,45]
elif region == "JS":
    box = [130,150,30,60]
else:
    raise ValueError("Unknown region")

date = sys.argv[1]
# date = '2099-01-17-18'
# year = '2099'
output_dir = sys.argv[2]
# output_dir = f'/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/lagranto/cao_input/{region}/{year}'
fig_path = '/home/nfieldhouse/ma/plots/'
format_in = "%Y-%m-%d-%H"
datetime = dt.datetime.strptime(date, format_in)
year = datetime.year
print(date)
print(year)
print(f"Region: {region}, Output Directory: {output_dir}")

# %%
# -----------------------------------
# load data from CESM
# -----------------------------------
timeframe = sys.argv[4]
# timeframe = 'future'
if timeframe == "present":
    cesm_type = 'HIST'
elif timeframe == "future":
    cesm_type = 'SSP370'
cesm = f"/net/meso/climphys/cesm212/b.e212.B{cesm_type}cmip6.f09_g17.{member_id}/archive/atm/hist/"
h3 = xr.open_dataset(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.{member_id}.cam.h3.{year}-01-01-00000.nc") # primary variables

# get dimensions
lons = h3.lon
lats = h3.lat
time = h3.time

# extract data from target date
h3["time"] = h3.indexes["time"].to_datetimeindex() # set index time of the dataarray to datetime
h3_date = h3.sel(time=datetime)

# get variables
PS = h3_date.variables["PS"].data  # surface pressure 
T = h3_date.variables["T"].data[25,:,:]  # temperature (3D on several pressure levels) on the 25th level, corresponds to plevel=859hPa
SST = h3_date.variables["SST"].data # sea surface temperature 

# calculate potential temperature at 850hPa and sea surface potential temperature (skintemperature)
SSTH = SST*(1000/(PS/100))**(0.286)
TH_850 = T*(1000/(85953/100))**(0.286)

# calculate 0,1-masks
CAOI = (SSTH - TH_850)
CAO_mask = np.where(CAOI > 8, 1,0) # threshold level 8K from casestudy

# exclude grid points with sea ice  (SST < 273K) & land
sea_ice = np.where(SST < 273, 0,1)
LSM = np.where(SST < 271,0,1)

# select only masks that are within our region
lonmin = box[0]
lonmax = box[1]
latmin = box[2]
latmax = box[3]
grid = np.meshgrid(lons,lats)
lon_mask_min = np.where(grid[0] > lonmin, 1, 0)
lon_mask_max = np.where(grid[0] < lonmax, 1, 0)
lat_mask_min = np.where((grid[1] > latmin), 1, 0)
lat_mask_max = np.where((grid[1] < latmax), 1, 0)
box_mask = lat_mask_min*lat_mask_max*lon_mask_min*lon_mask_max
mask = CAO_mask*box_mask*LSM*sea_ice

# %%
# save mask as a 0,1-field
fout = f"{output_dir}/cao_{date}"
field_out = xr.DataArray(name = "mask", data = mask, dims = ["lat", "lon"])
data_var = {"mask" : field_out}
ds_new = xr.Dataset(coords = {"lon": lons, "lat": lats}, data_vars = data_var)
ds_new.to_netcdf(fout)

# %%
# ----------------------------------------------
# make startf with upper boundary being PBL (for Lagranto)
# ----------------------------------------------
# assume PBL height of 850hPa (this is what we have roughly in cases in ERA5), note: we don't have PBL height on a high enough resolution
lower_boundary = 1000
upper_boundary = 850

n_levels = int((lower_boundary-upper_boundary)/25) + 1

# -------------------------- get list of gridpoints as input into LAGRANTO with PBL height ------------------------------------------
gp_lagranto = np.where(mask == 1)
gp_lat = gp_lagranto[0]*0.9424080000000004-90  # output of np.where are indices (not lons directly)
gp_lon = gp_lagranto[1]*1.25

# repeat lons and lats for number of height levels
lon_rep = np.repeat(gp_lon,n_levels)
lat_rep = np.repeat(gp_lat,n_levels)

p_rep = np.empty(0)
for i in range(len(gp_lon)):
    p_arr = np.arange(lower_boundary, upper_boundary-25, -25) # make array every 25hPa from lower to upper boundary
    p_rep = np.concatenate((p_rep, p_arr), axis=0)

time_rep = np.zeros(len(lon_rep))

df = pd.DataFrame(list(zip(time_rep, lon_rep, lat_rep, p_rep)), columns = ["time", "lon","lat","plev"])
fout = f"{output_dir}/CAO_{date}"
df.to_csv(fout, sep = " ", index = False, header = False)
# %%

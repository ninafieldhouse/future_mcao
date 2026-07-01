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
from matplotlib.patches import Rectangle
import cmcrameri.cm as cmc



# @jit(nopython=True, boundscheck=False, cache=True, fastmath=True)
# def v_interpolate_loop(z_out, data_out, z_in, data_in):
#     ny, nx = z_out.shape
#     for i in range(nx):
#         for j in range(ny):
#             data_out[j, i] = np.interp(
#                 z_out[j, i],
#                 z_in[:, j, i],
#                 data_in[:, j, i],
#             )

# def vertical_interpolate_3d(z_out, z_in, data_in):
#     ny, nx = z_out.shape
#     data_out = np.zeros((ny, nx))
#     v_interpolate_loop(z_out, data_out, z_in, data_in)
#     return data_out

# Define a function to create a custom colormap with a white segment
def create_custom_colormap(cmap, vmin, vmax, white_range):
    colors = cmap(np.linspace(0, 1, 256))
    white = np.array([1, 1, 1, 1])
    norm = plt.Normalize(vmin, vmax)
    white_min, white_max = norm(white_range[0]), norm(white_range[1])
    
    # Create new color array with white segment
    new_colors = np.zeros((256, 4))
    for i, val in enumerate(np.linspace(0, 1, 256)):
        if val < white_min:
            new_colors[i] = colors[int(val / white_min * 128)]
        elif val > white_max:
            new_colors[i] = colors[min(int((val - white_max) / (1 - white_max) * 128) + 128, 255)]  # Ensure index does not exceed 255
        else:
            new_colors[i] = white
    
    return LinearSegmentedColormap.from_list('custom_cmap', new_colors)

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

# %%
# -----------------------------------
# define boundary conditions
# -----------------------------------
# region = sys.argv[3]
region = 'JS'
if region == "NA":
    extent_map = [275,310,20,50] # for map
    box = [285, 295, 30, 45] # for mask
elif region == "JS":
    extent_map = [115,180,20,60] # for map
    box = [130,150,30,60] # for mask
else:
    raise ValueError("Unknown region")

box = {
    "JS": [130, 150, 30, 60],  # [lon_min, lon_max, lat_min, lat_max]
    "NA": [285, 295, 30, 45]   # [lon_min, lon_max, lat_min, lat_max]
}

# date = sys.argv[1]
date = '2091-02-22-12' # 1000/JS
member_id = 1400
# output_dir = sys.argv[2]
output_dir = f'/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/lagranto/cao_input/{region}'
data_dir = f'/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/cao/clim/{region}'
fig_path = '/home/nfieldhouse/ma/plots/casestudy'
format_in = "%Y-%m-%d-%H"
datetime = dt.datetime.strptime(date, format_in)
year = datetime.year
print(date)
print(year)

# %%
# -----------------------------------
# load data from CESM
# -----------------------------------
# timeframe = sys.argv[4]
timeframe = 'future'
if timeframe == "present":
    cesm_type = 'HIST'
elif timeframe == "future":
    cesm_type = 'SSP370'
cesm = f"/net/meso/climphys/cesm212/b.e212.B{cesm_type}cmip6.f09_g17.{member_id}/archive/atm/hist/"
h3 = xr.open_dataset(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.{member_id}.cam.h3.{year}-01-01-00000.nc") # primary variables
h4 = xr.open_dataset(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.{member_id}.cam.h4.{year}-01-01-00000.nc") # secondary variables


# get dimensions
lons = h3.lon
lats = h3.lat
time = h3.time

# extract data from target date
h3["time"] = h3.indexes["time"].to_datetimeindex() # set index time of the dataarray to datetime
h3_date = h3.sel(time=datetime, lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3]))
h4["time"] = h4.indexes["time"].to_datetimeindex() # set index time of the dataarray to datetime
h4_date = h4.sel(time=datetime, lon=slice(box[region][0], box[region][1]), lat=slice(box[region][2], box[region][3]))

# get variables
SST = h3_date.SST
print('Opened SST')
T = h3_date.T[25,:,:]
print('Opened T')
PS = h3_date.PS
print('Opened PS')
SHFLX = h3_date.SHFLX
print('Opened SHFLX')
LHFLX = h3_date.LHFLX
print('Opened LHFLX')


actual_sea_ice =  h4_date.ICEFRAC

SSTH = SST*(1000/(PS/100))**(0.286)
TH_850 = T*(1000/(85953/100))**(0.286)

# exclude grid points with sea ice  (SST < 273K) & land
sea_ice_used_in_stats = np.where(SST < 273, 0,1)
LSM = np.where(SST < 271,0,1)
sea_ice_approx = LSM+ sea_ice_used_in_stats # approximate sea ice, used for CAOI calculation
mask = np.where(sea_ice_approx== 1, 1,0) # mask for land and sea ice
CAOI = SSTH - TH_850
# calculate turbulent heat flux
THFLX = SHFLX + LHFLX

sea_ice = np.where(SST < 273, np.nan,1)
LSM = np.where(SST < 271,np.nan,1)
mask = LSM * sea_ice

# calculate 0,1-masks
CAOI = (SSTH - TH_850)* mask

CAO_mask = np.where(CAOI > 4, 1,0)
# %%
# Plot bias between mask and actual sea ice
lons = h3_date.lon
lats = h3_date.lat
bias = mask - actual_sea_ice

cmap_name = 'davos_r'
boundaries = np.arange(0, 19, 1)
cbarticks = np.arange(0, 19, 2)

cmap, n_main = cmap_pcolormesh(cmap_name, boundaries) #main colors, without extended colors
norm = BoundaryNorm(boundaries, ncolors=n_main, clip=False)


fig, ax = plt.subplots(1, 1, figsize=(7, 5), subplot_kw={'projection': ccrs.PlateCarree()})

# pcm = ax.contourf(lons, lats, CAOI, transform=ccrs.PlateCarree(), cmap=cmap, levels=boundaries, extend='both')
im = ax.pcolormesh(lons, lats, CAOI, shading='auto',
                       transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
ax.coastlines()
ax.set_title(f'Bias: Mask - Actual Sea Ice\nEnsemble: {member_id}, Date: {date}', fontsize=14)
ax.set_extent([box[region][0], box[region][1], box[region][2], box[region][3]], crs=ccrs.PlateCarree())

cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.05, extend ='both')
cbar.set_label('Bias (Mask - Actual Sea Ice)')

plt.savefig(f"{fig_path}/CAOI_{region}_{date}.png", bbox_inches='tight', dpi=200, facecolor='white', transparent=False)
plt.show()
# %%



#%%
# -----------------------------------
# case study threshold levels
# -----------------------------------
# value = float(sys.argv[5])
value = 5
if value >= 5:
    thresholds = np.array([4, 6, 8, 12])
    # thresholds = np.array([4])
    fig, axs = plt.subplots(3, len(thresholds), 
                            figsize=(15, 10),
                            subplot_kw={'projection': ccrs.PlateCarree()}, layout='constrained')

    for i, threshold in enumerate(thresholds):
        # Plot CAOI with hatches in the first row
        CAO_mask = np.where(CAOI > threshold, 1, 0)
        mask = CAO_mask * box_mask * LSM * sea_ice
        extended_lonmin = lonmin 
        extended_lonmax = lonmax 
        extended_latmin = latmin 
        extended_latmax = latmax 

        # Create a copy of CAOI for plotting and set all values <= 0 to NaN
        CAOI_plot = CAOI.copy()
        CAOI_plot[CAOI_plot <= 0] = np.nan

        # Use the colormap, levels, and colorbar settings from CAO_casestudy.py
        boundaries = np.arange(0, 21, 2)  # Match levels from CAO_casestudy.py
        norm = BoundaryNorm(boundaries, ncolors=256)
        cmap = plt.get_cmap('plasma')  # Match colormap from CAO_casestudy.py

        # Plot CAOI
        ax = axs[0, i]
        ax.coastlines()
        ax.set_extent([extent_map[0], extent_map[1], extent_map[2], extent_map[3]], crs=ccrs.PlateCarree())
        contour_CAOI = ax.contourf(lons, lats, CAOI_plot, transform=ccrs.PlateCarree(), cmap=cmap, levels=boundaries, extend='max')

        # # Add hatching for areas where mask > 0
        # hatch_mask = np.where(mask > 0, 1, np.nan)
        # ax.contourf(lons, lats, hatch_mask, transform=ccrs.PlateCarree(), colors='none', hatches=['/////'], levels=[0, 0.5, 1], extend='neither')
        
        lon_min, lon_max = 130, 150
        lat_min, lat_max = 30, 60

        # Create a rectangle with the specified coordinates
        box = Rectangle(
            (lon_min, lat_min),  # Bottom-left corner (longitude, latitude)
            lon_max - lon_min,   # Width (longitude range)
            lat_max - lat_min,   # Height (latitude range)
            linewidth=2,         # Line width of the box
            edgecolor='black',      # Color of the box edges
            facecolor='none',     # Transparent fill
            transform=ccrs.PlateCarree(),  # Ensure the box is in the correct projection
            zorder=3              # Draw on top of the trajectories
        )

        # Add the rectangle to the first plot
        ax = axs[0,0]
        ax.add_patch(box)

        # Plot turbulent heat flux in the second row
        cmap = plt.get_cmap('coolwarm')  # Keep the custom colormap for heat flux
        vmin, vmax = -1500, 1500
        white_range = (-300, 300)
        custom_cmap = create_custom_colormap(cmap, vmin, vmax, white_range)
        boundaries = np.arange(vmin, vmax + white_range[1], white_range[1])
        norm = BoundaryNorm(boundaries, ncolors=256)

        ax = axs[1, i]
        ax.coastlines()
        ax.set_extent([extent_map[0], extent_map[1], extent_map[2], extent_map[3]], crs=ccrs.PlateCarree())
        contour_THFLX = ax.contourf(lons, lats, THFLX, transform=ccrs.PlateCarree(), cmap=custom_cmap, levels=boundaries, extend='both')

        # Plot contour lines of mask with different threshold levels
        for t in thresholds:
            linewidth = 1.5 if t == threshold else 0.5
            ax.contour(lons, lats, CAOI, levels=[t], colors='black', linewidths=linewidth, transform=ccrs.PlateCarree())

        # Plot CAO mask in the third row
        mask_cmap = ListedColormap(['white', 'cornflowerblue'])  # Binary colormap: white for 0, blue for 1

        ax = axs[2, i]
        ax.coastlines()
        ax.set_extent([extent_map[0], extent_map[1], extent_map[2], extent_map[3]], crs=ccrs.PlateCarree())
        ax.contourf(lons, lats, mask, transform=ccrs.PlateCarree(), cmap=mask_cmap, levels=[0, 0.5, 1], extend='neither')

        # Add threshold level to the column title
        axs[0, i].set_title(f'Threshold = {threshold}')

    # Add main title
    fig.suptitle(f'CAOI and CAO Mask on {date} for region {region}', fontsize=16)

    # Add colorbar for CAOI
    cbar = fig.colorbar(contour_CAOI, ax=axs[0, :], orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('CAOI')
    cbar.set_ticks(np.arange(0, 21, 2))  # Match ticks from CAO_casestudy.py

    # Add colorbar for turbulent heat flux
    cbar = fig.colorbar(contour_THFLX, ax=axs[1, :], orientation='horizontal', pad=0.05, aspect=50)
    cbar.set_label('Turbulent Heat Flux')
    cbar.set_ticks(boundaries)
    # plt.show()
    plt.savefig(f"{fig_path}/casestudy/CAOI_and_Mask_{region}_{date}.png", bbox_inches='tight', dpi=200, facecolor='white', transparent=False)
# %%

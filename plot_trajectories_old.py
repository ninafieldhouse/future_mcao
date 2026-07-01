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


def findLongestConseqSubseq_time(arr, n):
    #The longest sequence in an empty array is, of course, 0, so we can simply return that. 
    if n==0:
       return 0
 
    # check each possible sequence from the start
    # then update optimal length
    ans = 1
    count = 1
    t= []
    for i in range(1, n):
      #For handling duplicate elements
      if arr[i]!=arr[i-1]:
 
        # if difference between current element and previous element is 1
        # then we want to update our current sequence count
 
        if arr[i]-arr[i-1] == 1:
            count += 1
            t.append(i)
 
        # otherwise, we will update our count to zero to check for other sequences.
        # before updating count value we have to check if current sequence length is more than our ans.
        # if count > ans then we want to update our ans.
 
        else:
            ans = max(ans, count)
            count = 1

    # print(t[-1]+1)
 
    return max(ans,count), t[-1]+1 #To handle the case in which last element is present in longest sequence.

def findLongestConseqSubseq(arr, n):
    #The longest sequence in an empty array is, of course, 0, so we can simply return that. 
    if n==0:
       return 0
 
    # check each possible sequence from the start
    # then update optimal length
    ans = 1
    count = 1

    for i in range(1, n):
      #For handling duplicate elements
      if arr[i]!=arr[i-1]:
 
        # if difference between current element and previous element is 1
        # then we want to update our current sequence count
 
        if arr[i]-arr[i-1] == 1:
            count += 1
 
        # otherwise, we will update our count to zero to check for other sequences.
        # before updating count value we have to check if current sequence length is more than our ans.
        # if count > ans then we want to update our ans.
 
        else:
            ans = max(ans, count)
            count = 1
 
    return max(ans,count) #To handle the case in which last element is present in longest sequence.

# %%
# -------------------------------
#  plot backward trajectories
# -------------------------------

file = "B_19901207_12_with_TH.txt"
path = "/net/litho/atmosdyn2/nfieldhouse/data/1000/lagranto/traj/JS/"
all_traj = pd.read_csv(path+file)
trajs = Tra(path+file)
 
# # -------------------------
# plotvar = 'p'    # color trajetories based on their pressure
# unit    = '[hPa]' # important for plotting
# cmap = plt.get_cmap('Spectral')
# levels  = np.arange(400,1051,50)
# # -------------------------
# -------------------------
plotvar = 'TH'    # color trajetories based on their pressure
unit    = '[K]' # important for plotting
cmap = plt.get_cmap('Spectral_r')
levels  = np.arange(230,301,5)
# -------------------------

# # -------------------------
# plotvar = 'SST'    # color trajetories based on their pressure
# unit    = '[K]' # important for plotting
# cmap = plt.get_cmap('Spectral_r')
# levels = np.arange (260, 301, 5)
# # -------------------------

slp_levels = np.arange(958, 1070, 4)
map_projection = ccrs.PlateCarree(central_longitude=180)
map_projection = ccrs.AzimuthalEquidistant(central_longitude=120, central_latitude=55)

# %%
fig, ax = plt.subplots(1, 1, figsize = (15,12), subplot_kw = dict(projection = map_projection), layout = 'constrained') #create subplots
ax.add_feature(cfeature.LAND, facecolor = (217/256, 217/256, 217/256))
ax.set_extent([60,260,25,90], crs = ccrs.PlateCarree()) # extent as box: [115,180,20,60]
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
# gl.bottom_labels = False
gl.xlabel_style = {"size" : 20, "color" : "k"}
gl.ylabel_style = {"size" : 20, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

# Define the coordinates of the box
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
# Add the rectangle to the plot
ax.add_patch(box)

k = plot_trajs(ax,trajs[::-5,:],plotvar,cmap=cmap,alpha = 0.7,levels=levels,linewidth=2) # trajs[::-10,:] plots every 10th trajectory starting from top, all timesteps of the trajectory
# k = plot_trajs(ax,trajs[::-50,:],plotvar,cmap=cmap,alpha = 0.7,levels=levels,linewidth=2) # plot the first trajectory

# Adjust the colorbar position to avoid overlap
# cax = fig.add_axes([0.15, 0.05, 0.7, 0.03])  # Move the colorbar below the plot
cbar = plt.colorbar(k, ax=ax, orientation="horizontal", shrink =0.6, pad=0.05)
cbar.set_label(label=f'{plotvar} {unit}', size=20)
cbar.ax.tick_params(labelsize=15)

# ax.clabel(k, fontsize = 10) #replaced m with k

# Save the figure
fig_path = '/home/nfieldhouse/ma/plots/'
fig.savefig(f'{fig_path}traj/{file}_{plotvar}_test.png', bbox_inches = "tight")
plt.close(fig)

# %%
# -------------------------------
# plot SST trajectories over sea ice
# -------------------------------

levels = np.array([0, 265, 270, 273, 278, 283, 288, 293, 298, 303])
timeframe = 'present'
year = '1997'
date = '1997-01-25-06'
format_in = "%Y-%m-%d-%H"
datetime = dt.datetime.strptime(date, format_in)

if timeframe == "present":
    cesm_type = 'HIST'
elif timeframe == "future":
    cesm_type = 'SSP370'
cesm = f"/net/meso/climphys/cesm212/b.e212.B{cesm_type}cmip6.f09_g17.1000/archive/atm/hist/"
h3 = xr.open_dataset(cesm + f"b.e212.B{cesm_type}cmip6.f09_g17.1000.cam.h3.{year}-01-01-00000.nc") # primary variables
# get dimensions
lons = h3.lon
lats = h3.lat
time = h3.time

# extract data from target date
h3["time"] = h3.indexes["time"].to_datetimeindex() # set index time of the dataarray to datetime
h3_date = h3.sel(time=datetime)
SST = h3_date.variables["SST"].data # sea surface temperature

# Get the longitude and latitude grid
lon_grid, lat_grid = np.meshgrid(lons, lats)

# Get the land geometries from Cartopy's Natural Earth dataset
land_feature = cfeature.NaturalEarthFeature('physical', 'land', '110m')
land_geometries = list(land_feature.geometries())

# Prepare the geometries for faster computation
prepared_land = [prep(geom) for geom in land_geometries]

# Initialize the land mask (1 for water, 0 for land)
sea_mask = np.ones(lon_grid.shape, dtype=int)

lon_grid = np.where(lon_grid > 180, lon_grid - 360, lon_grid)

# Loop through each grid point and check if it falls on land
for i in range(lon_grid.shape[0]):
    for j in range(lon_grid.shape[1]):
        point = Point(lon_grid[i, j], lat_grid[i, j])
        if any(prepared_geom.contains(point) for prepared_geom in prepared_land):
            sea_mask[i, j] =0  # Mark as land

# Now `land_mask` contains 0 for land and 1 for water

    
# Create a sea ice mask (1 where SST < 273, 0 otherwise)
threshold = 273  # Kelvin
ice_mask = np.where(SST < threshold, 1, 0)
sea_ice_mask = ice_mask * sea_mask  # Combine with land mask

# Initialize a list to store indices of trajectories that pass over sea ice
sea_ice_traj_indices = []

# Loop through each trajectory
for traj_idx in range(len(trajs)):
    traj_lons = trajs["lon"][traj_idx, :]  # Longitudes of the trajectory
    traj_lats = trajs["lat"][traj_idx, :]  # Latitudes of the trajectory

    # Check if any point in the trajectory passes over sea ice
    for lon, lat in zip(traj_lons, traj_lats):
        # Find the nearest grid point
        lon_idx = np.abs(lons - lon).argmin()
        lat_idx = np.abs(lats - lat).argmin()

        # Check if the grid point is over sea ice and not over land
        if sea_ice_mask[lat_idx, lon_idx] == 1:
            sea_ice_traj_indices.append(traj_idx)
            break  # No need to check further points for this trajectory

# Filter the trajectories to include only those that pass over sea ice
sea_ice_trajs = Tra()
sea_ice_trajs.set_array(trajs[sea_ice_traj_indices, :])

# Plot the trajectories that pass over sea ice
fig, ax = plt.subplots(1, 1, figsize=(15, 12), subplot_kw=dict(projection=map_projection), layout = 'constrained')
ax.set_extent([60, 260, 25, 90], crs=ccrs.PlateCarree())
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels=True, dms=True)
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {"size": 20, "color": "k"}
gl.ylabel_style = {"size": 20, "color": "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
gl.ylocator = mticker.FixedLocator(np.arange(0, 100, 20))


sea_ice_cmap = ListedColormap(['white', 'cornflowerblue'])  # White for no sea ice, blue for sea ice
sea_ice_levels = [0, 0.5, 1]  # Binary levels for the mask
sea_ice_plot = ax.contourf(lons, lats, sea_ice_mask, levels=sea_ice_levels, cmap=sea_ice_cmap, transform=ccrs.PlateCarree(), zorder=1)
ax.add_feature(cfeature.LAND, facecolor=(217/256, 217/256, 217/256), zorder =2)

# Plot the filtered trajectories
k = plot_trajs(ax, sea_ice_trajs[::-10, :], plotvar, cmap=cmap, alpha=0.7, levels=levels, linewidth=2)

# Define the coordinates of the box
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

# Add the rectangle to the plot
ax.add_patch(box)

# Add a colorbar
cbar = plt.colorbar(k, ax=ax, orientation="horizontal", shrink=0.6, pad=0.05)
cbar.set_label(label=f'{plotvar} {unit}', size=20)
cbar.ax.tick_params(labelsize=15)

# Save the figure
fig.savefig(f'{fig_path}traj/{file}_sea_ice_{threshold}_{plotvar}.png', bbox_inches="tight", dpi = 200)
plt.close(fig)
# %%
# -------------------------------
#  identify different air streams
# -------------------------------
###### select BCAR trajectories
# select that p always larger than 700hPa
low_inds = np.where(trajs["p"]>700)
inds = []
for j in range(len(trajs)):
    count  = len(np.where(low_inds[0] == j)[0])     # get those indices where condition is met for all timesteps
    if count == 29:
        inds.append(j)
low_trajs = Tra()
low_trajs.set_array(trajs[inds,:])

# create array with indices for low level trajectoris BUT non BCAR trajetories
array = np.arange(0, len(low_trajs),1)

# select for temperature threshold during 3 days
cold_inds = np.where((low_trajs["TH"]<280) & (low_trajs["SST"]==0))
inds_traj = np.unique(cold_inds[0])
inds = []
for ind in inds_traj:
    traj_i = np.where(cold_inds[0] == ind)[0]
    timesteps = []
    for j in traj_i:
        time = cold_inds[1][j]
        timesteps.append(time)
    consec = findLongestConseqSubseq(timesteps, len(timesteps))
    if consec > 12:
        inds.append(ind)

cold_trajs = Tra()
cold_trajs.set_array(low_trajs[inds,:])

#  remove non cold trajs from array for non BCAR trajetories
if len(cold_trajs) != len(low_trajs):
    non_cold = np.delete(array, inds)
else: 
    non_cold = array

# select that at least at one timestep p has to be higher than 1020hPa
anti_inds = np.where(cold_trajs["PSL"]>1020)
anti_inds = np.unique(anti_inds[0])
anti_trajs = Tra()
anti_trajs.set_array(cold_trajs[anti_inds,:])

if len(cold_trajs) != len(anti_trajs):
    print("no anti inds " + date)

non_BCAR_trajs = Tra()
non_BCAR_trajs.set_array(low_trajs[non_cold,:])

BCAR_trajs = anti_trajs

# fraction of low level trajectories
ll_inds = np.where(trajs["p"]>700)
inds = []
for j in range(len(trajs)):
    count  = len(np.where(ll_inds[0] == j)[0])     # get those indices where condition is met for all timesteps
    if count == 169:
        inds.append(j)
ll_trajs = Tra()
ll_trajs.set_array(trajs[inds,:])

# if len(BCAR_trajs) != 0:
#     BCAR_fract = len(anti_trajs)/len(ll_trajs)
#     BCAR_fract_overall = len(anti_trajs)/len(trajs)
#     ll_fract_overall = len(ll_trajs)/len(trajs)

# else:
#     BCAR_fract = 0
#     BCAR_fract_overall = 0
#     ll_fract_overall = len(ll_trajs)/len(trajs)

###### select trajectories that make upper level contribution
ll_inds = np.where(trajs["p"]>700)
inds = []
for j in range(len(trajs)):
    count  = len(np.where(ll_inds[0] == j)[0])     # get those indices where condition is met for all timesteps
    if count == 169:
        inds.append(j)
ll_trajs = Tra()
ll_trajs.set_array(trajs[inds,:])

all = np.where(trajs[:,0])[0]
mask = np.isin(all, inds, invert = True)
ul_inds = all[mask]

ul_trajs = Tra()
ul_trajs.set_array(trajs[ul_inds,:])

DI = np.where(abs((trajs["p"][:,:1] - trajs["p"][:,:12]))>400)
DI_unique = np.unique(DI[0])

DI_trajs = Tra()
DI_trajs.set_array(trajs[DI_unique,:])

# if len(DI_trajs) != 0:
#     DI_fract = len(DI_trajs)/len(ul_trajs)
#     DI_fract_overall = len(DI_trajs)/len(trajs)
#     ul_fract_overall = len(ul_trajs)/len(trajs)

# else:
#     DI_fract = 0
#     DI_fract_overall = 0
#     ul_fract_overall = len(ul_trajs)/len(trajs)

weird = np.unique(np.where(DI_trajs["TH"]<0)[0])
DI = np.where(DI_trajs[:,0])[0]
good_inds = np.isin(DI, weird, invert = True)
good = DI[good_inds]

DI_tras_good = Tra()
DI_tras_good.set_array(DI_trajs[good,:])

fout_BCAR = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_BCAR_trajs_" + date
BCAR_trajs.write_ascii(fout_BCAR)

if len(DI_tras_good) != 0:
    fout_DI = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_DI_trajs_" + date
    DI_tras_good.write_ascii(fout_DI)

fout_non_BCAR = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_nonBCAR_low_level_trajs_" + date
non_BCAR_trajs.write_ascii(fout_non_BCAR)

fout_ll = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_all_low_level_trajs_" + date
ll_trajs.write_ascii(fout_ll)
# %%

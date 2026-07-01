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
import sys
from geopy.distance import great_circle
from dypy.small_tools  import interpolate
import xarray as xr
import random
from datetime import datetime, timedelta
from tqdm import tqdm
from multiprocessing import Pool
import glob
from scipy import stats
from cartopy.util import add_cyclic_point

plt.style.use('default') # I personally prefer seaborn for the graph style, but you may choose whichever you want.
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "grid.color" :  (217/256, 217/256, 217/256), #(217/256, 217/256, 217/256), "slategrey"
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)

def generate_random_dates(start_year, end_year, num_dates, rng):

    # Initialize an empty list to store the dates
    valid_dates = []

    # Generate all potential valid dates within the range
    potential_dates = []
    for year in range(start_year, end_year + 1):
        for month in [12, 1, 2]:  # December, January, February
            if month == 1 or month == 2:
                year_offset = year + 1 if month in [1, 2] else year
            else:
                year_offset = year
            if month == 12:
                next_month_date = datetime(year_offset + 1, 1, 1)  # Move to January of the next year
            else:
                next_month_date = datetime(year_offset, month + 1, 1)
            days_in_month = (next_month_date - timedelta(days=1)).day
            for day in range(1, days_in_month + 1):
                for hour in range(24): 
                    potential_dates.append(datetime(year_offset, month, day, hour))

    # Randomly select dates ensuring 5 days apart
    while len(valid_dates) < num_dates:
        random_date = rng.choice(potential_dates)
        if all(abs((random_date - date).days) >= 5 for date in valid_dates):
            valid_dates.append(random_date)

    # Sort dates before returning
    valid_dates.sort()
    return valid_dates

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
   pvalues=np.array(pvalues)
   pthresh = fdr(pvalues.reshape(-1),alpha)
   return pvalues<pthresh

# %%
var = "jet"
time = "96"
for region in ["Siberia"]:
        print(region)

        for i in range(3):  # loop through groups 
            print(i)
            
            if i < 2:  
                folder_path = f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/significance/{region}/{var}"   # Replace with your folder path
                file_pattern = f"comp_anom_{var}_group_{i+1}_CAO_{time}hafter_{region}_*"  # Adjust this pattern if needed (e.g., based on file names)
            elif i == 2:
                folder_path = f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/significance/{region}/{var}"   # Replace with your folder path
                file_pattern = f"comp_anom_{var}_group_1_CAO_{time}hafter_{region}_*"  # Adjust this pattern if needed (e.g., based on file names)

            # Get a list of all NetCDF files
            file_list = sorted(glob.glob(f"{folder_path}/{file_pattern}"))  # Sort for consistent order

            # Load all NetCDF files into an xarray dataset
            datasets = [xr.open_dataset(file) for file in file_list]

            # Stack the 2D fields into a single 3D array (time, lat, lon)
            # Assume the 2D field in each file has the same variable name, e.g., 'data'
            data_array = xr.concat([ds['anomaly'] for ds in datasets], dim='time')

            # read in composite of anomalies
            if region == "Siberia":
                ds_var = xr.open_dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/composites/comp_anom_{var}_group_{i+1}_CAO_{time}hafter_{region}") # SLP
            
            anom = ds_var.variables["anomaly"]
            
            # find quantile/p-value of anomalies
            percentiles = np.zeros(anom.shape)  # To store results

            # Loop through each grid point
            for y in range(anom.shape[0]):  # Loop over rows 
                for z in range(anom.shape[1]):  # Loop over columns
                    score = anom[y, z]
                    dist = data_array[:, y, z]  # Distribution at this grid point
                    
                    # Compute percentile
                    percentiles[y, z] = np.mean(dist <= score)

            # those at right tale (1-percentile) to make them one-sided
            percentiles_one_sided = np.where(percentiles > 0.5, 1-percentiles, percentiles)

            # call FDR function
            fdr_out = fdr_xr_2d(pvalues = percentiles_one_sided, alpha = 0.1)
            masked = xr.DataArray(name = "masked",data = fdr_out * anom, dims=anom.dims)

            fout = f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/composites/comp_anom_{var}_group_{i}_CAO_{time}hafter_{region}_significant_fdr_alpha_0.1"
            data_var = {"masked": masked}
            ds_new = xr.Dataset(coords = ds_var.coords, data_vars = data_var)
            ds_new.to_netcdf(fout)                
            
# %%
# plot_var = "all_cyclone_freq_anom"
# # region = "Siberia"
# region = "NA"

# # titles = ["CAO - 8days","CAO - 6days","CAO - 4days", "CAO  - 2days", "CAO"]
# row_headers = ["small BCAR contribution", "mean BCAR contribution", "large BCAR contribution"]

# # cmap_th
# cmap = matplotlib.cm.get_cmap("binary")
# white = cmap(0)
# cmap = matplotlib.cm.get_cmap("RdBu_r")
# cmap_start = []
# for i in np.linspace(0,0.45,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)
# cmap_start.append(white)
# for i in np.linspace(0.55,1.01,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)    
# cmap_colors = cmap_start
# cmap_th = ListedColormap(cmap_colors, 'cmap_temp')


# # plot composites
# if region == "Siberia":
#     extent = [120,220,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude=180)
#     xy = (-185,40)

# elif region == "NA":
#     extent = [-160,20,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude= 0)
#     xy = (-170,40)

# fig = plt.figure(figsize = (50,15))
# gs = gridspec.GridSpec(1,3, wspace = 0.08)

# for i in range(3):                  # plot each group

#         # read in data
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_{region}_group_{i+1}_0_to_3days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         var = ds_var.variables["masked"][:].data
#         lons = np.array(ds_var.variables["lon"])
#         lats = np.array(ds_var.variables["lat"])

#         ax = fig.add_subplot(gs[0,i], projection = map_projection)

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")

#         if i == 0:
#             ax.set_title(row_headers[i], fontsize = 30)

#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "black", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         gl.bottom_labels = False
#         gl.xlabel_style = {"size" : 20, "color" : "k"}
#         gl.ylabel_style = {"size" : 20, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
    
#         h = ax.contourf(np.arange(-180,180,0.5), np.arange(-90,90.5,0.5), var[:,:], cmap = cmap, levels = np.arange(-0.07,0.08,0.02), extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)           
        

# cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.02])
# cbar = fig.colorbar(h, cax=cbar_ax, orientation= "horizontal")
# cbar.set_label(label = "frequency anomaly [%]", size = 30)
# cbar.ax.tick_params(labelsize = 25)

#  %%

# plot_var = "freq_anom"
# region = "Siberia"
# # region = "NA_sea_ice_allowed"

# titles = ["small BCAR contribution", "mean BCAR contribution", "large BCAR contribution"]

# # cmap_sst
# cmap = matplotlib.cm.get_cmap("binary")
# white = cmap(0)
# cmap = matplotlib.cm.get_cmap("RdBu_r")
# cmap_start = []
# for i in np.linspace(0,0.45,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)
# cmap_start.append(white)
# for i in np.linspace(0.55,1.01,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)    
# cmap_colors = cmap_start
# cmap_sst = ListedColormap(cmap_colors, 'cmap_temp')

# # plot composites
# if region == "Siberia":
#     extent = [120,240,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude=180)
#     xy = (80,40)

# elif region == "NA_sea_ice_allowed":
#     extent = [-120,20,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude= 0)
#     xy = (-90,40)

# cmap_anom = cmap_sst
# levels_anom = np.arange(-0.09,0.11,0.02)
# norm_anom = matplotlib.colors.TwoSlopeNorm(vmin = -0.1, vmax = 0.1, vcenter = 0)
# color_label_anom = r"cyclone frequency anomaly"

# # cmap_sst
# cmap = matplotlib.cm.get_cmap("binary")
# white = cmap(0)
# cmap = matplotlib.cm.get_cmap("RdBu_r")
# cmap_start = []
# for i in np.linspace(0,0.45,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)
# cmap_start.append(white)
# for i in np.linspace(0.55,1.01,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)    
# cmap_colors = cmap_start
# cmap_sst = ListedColormap(cmap_colors, 'cmap_temp')

# # plot composites
# if region == "Siberia":
#     extent = [120,240,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude=180)
#     xy = (80,40)

# elif region == "NA_sea_ice_allowed":
#     extent = [-100,20,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude= 0)
#     xy = (-90,40)

# cmap_abs = "Oranges"
# levels_abs = np.arange(0,0.5,0.05)
# norm_abs = matplotlib.colors.TwoSlopeNorm(vmin = 0, vmax = 0.5, vcenter = 0.25)
# color_label_abs = r"cyclone frequency anomaly"


# path_clim = "/net/tropo/atmosdyn/mfederer/ape/analysis/clim/cmask_expl_mean_djf_1979_2019"
# clim_ds = nc.Dataset(path_clim, mode = "r")
# clim_expl = np.array(clim_ds.variables["INPUT_SEL"])[0,:,:]

# path_clim = "/net/tropo/atmosdyn/mfederer/ape/analysis/clim/cmask_mean_djf_1979_2019"
# clim_ds = nc.Dataset(path_clim, mode = "r")
# clim_all = np.array(clim_ds.variables["INPUT_SEL"])[0,:,:]

# fig = plt.figure(figsize = (27,12))
# gs = gridspec.GridSpec(2,3, wspace = 0.1)

# for k in range(3):              # plot each timestep

#     if k == 0:
#         # read in data
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anomaly_{region}_group{k+1}_3_to_5days", mode = "r") # SLP
#         var = np.array(ds_var.variables["anomaly"])
#         # if region == "Siberia":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_Siberia_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # elif region == "NA_sea_ice_allowed":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_NA_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # var_sign= ds_var_sign.variables["masked"][:].data
#         # lons = ds_var.variables["lon"][:].data
#         lons = np.arange(-180,180,0.5)
#         lats = np.arange(-90,90.5,0.5)
#         var, lons = add_cyclic_point(var, coord=lons)
#         # lats = ds_var.variables["lat"][:].data

#         # mask_neg = var_sign[:,:] < 0
#         # mask_pos = var_sign[:,:] > 0

#         # ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_{region}_group{k+1}_5days", mode = "r") # SLP
#         # var_cont = np.array(ds_var.variables["frequency"])
 
#         var = np.where(var >= 100, np.nan,var)

#         ax = fig.add_subplot(gs[0,k], projection = map_projection)
#         ax.set_title(titles[k], loc = "center", fontsize = 30)

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")

#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "gray", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         gl.bottom_labels = False
#         gl.xlabel_style = {"size" : 30, "color" : "k"}
#         gl.ylabel_style = {"size" : 30, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
#         # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
#         # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

#         h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =1)           
#         i = ax.contour(lons,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
#         # ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
#         # ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())         

#         # plot anomaly
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anomaly_{region}_group{k+1}_3to5days", mode = "r") # SLP
#         var = np.array(ds_var.variables["anomaly"])
#         # if region == "Siberia":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anom_Siberia_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # elif region == "NA_sea_ice_allowed":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anom_NA_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # var_sign= ds_var_sign.variables["masked"][:].data
#         # lons = ds_var.variables["lon"][:].data
#         lons = np.arange(-180,180,0.5)
#         lats = np.arange(-90,90.5,0.5)
#         var, lons = add_cyclic_point(var, coord=lons)
#         # lats = ds_var.variables["lat"][:].data

#         # mask_neg = var_sign[:,:] < 0
#         # mask_pos = var_sign[:,:] > 0

#         # ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_{region}_group{k+1}_0_to_3days", mode = "r") # SLP
#         # var_cont = np.array(ds_var.variables["frequency"])

#         var = np.where(var >= 100, np.nan,var)

#         ax = fig.add_subplot(gs[1,k], projection = map_projection)

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")

#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "gray", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         # gl.bottom_labels = False
#         gl.xlabel_style = {"size" : 30, "color" : "k"}
#         gl.ylabel_style = {"size" : 30, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
#         # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
#         # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

#         h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =1)           
#         i = ax.contour(lons,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
#         # ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
#         # ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())   

#     if k == 1:
#         # read in data
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anomaly_{region}_group{k+1}_3_to_5days", mode = "r") # SLP
#         var = np.array(ds_var.variables["anomaly"])
#         # if region == "Siberia":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_Siberia_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # elif region == "NA_sea_ice_allowed":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_NA_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # var_sign= ds_var_sign.variables["masked"][:].data
#         # lons = ds_var.variables["lon"][:].data
#         lons = np.arange(-180,180,0.5)
#         lats = np.arange(-90,90.5,0.5)
#         var, lons = add_cyclic_point(var, coord=lons)

#         # mask_neg = var_sign[:,:] < 0
#         # mask_pos = var_sign[:,:] > 0

#         # ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_{region}_group{k+1}_0_to_3days", mode = "r") # SLP
#         # var_cont = np.array(ds_var.variables["frequency"])

#         var = np.where(var >= 100, np.nan,var)

#         ax = fig.add_subplot(gs[0,k], projection = map_projection)

#         ax.set_title(titles[k], loc = "center",fontsize = 25)

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")


#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "gray", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         gl.bottom_labels = False
#         gl.left_labels = False
#         gl.xlabel_style = {"size" : 30, "color" : "k"}
#         gl.ylabel_style = {"size" : 30, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
#         # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
#         # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

#         h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =1)    
#         i = ax.contour(lons,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
#         # ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
#         # ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())   

#         # plot anomaly
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anomaly_{region}_group{k+1}_3to5days", mode = "r") # SLP
#         var = np.array(ds_var.variables["anomaly"])
#         # if region == "Siberia":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anom_Siberia_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # elif region == "NA_sea_ice_allowed":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anom_NA_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # var_sign= ds_var_sign.variables["masked"][:].data
#         # lons = ds_var.variables["lon"][:].data
#         # lats = ds_var.variables["lat"][:].data
#         lons = np.arange(-180,180,0.5)
#         lats = np.arange(-90,90.5,0.5)
#         var, lons = add_cyclic_point(var, coord=lons)

#         # mask_neg = var_sign[:,:] < 0
#         # mask_pos = var_sign[:,:] > 0
        
#         # ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_{region}_group{k+1}_0to3days", mode = "r") # SLP
#         # var_cont = np.array(ds_var.variables["frequency"])

#         var = np.where(var >= 100, np.nan,var)

#         ax = fig.add_subplot(gs[1,k], projection = map_projection)
        

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")

#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "gray", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         # gl.bottom_labels = False
#         gl.left_labels = False
#         gl.xlabel_style = {"size" : 30, "color" : "k"}
#         gl.ylabel_style = {"size" : 30, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
#         # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
#         # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

#         h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom  ,levels = levels_anom ,norm = norm_anom , extend = "both", transform=ccrs.PlateCarree(), alpha =1)    
#         i = ax.contour(lons,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
#         # ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
#         # ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())   

#     if k == 2:
#         # read in data
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anomaly_{region}_group{k+1}_3_to_5days", mode = "r") # SLP
#         var = np.array(ds_var.variables["anomaly"])
#         # if region == "Siberia":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_Siberia_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # elif region == "NA_sea_ice_allowed":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_anom_NA_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # var_sign= ds_var_sign.variables["masked"][:].data
#         # lons = ds_var.variables["lon"][:].data
#         lons = np.arange(-180,180,0.5)
#         lats = np.arange(-90,90.5,0.5)

#         # mask_neg = var_sign[:,:] < 0
#         # mask_pos = var_sign[:,:] > 0

#         # ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/cyclones/all_cyclone_freq_{region}_group{k+1}_0_to_3days", mode = "r") # SLP
#         # var_cont = np.array(ds_var.variables["frequency"])

#         var = np.where(var >= 100, np.nan,var)
#         var, lons = add_cyclic_point(var, coord=lons)

#         ax = fig.add_subplot(gs[0,k], projection = map_projection)
#         ax.set_title(titles[k], loc = "center", fontsize = 30)

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")

#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "gray", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         gl.bottom_labels = False
#         gl.left_labels = False
#         gl.xlabel_style = {"size" : 30, "color" : "k"}
#         gl.ylabel_style = {"size" : 30, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
#         # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
#         # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

#         h= ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =1)
#         i = ax.contour(lons,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
#         # ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
#         # ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())   


#         # plot anomaly
#         ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anomaly_{region}_group{k+1}_3to5days", mode = "r") # SLP
#         var = np.array(ds_var.variables["anomaly"])
#         # if region == "Siberia":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anom_Siberia_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # elif region == "NA_sea_ice_allowed":
#         #     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_anom_NA_group_{k+1}_5days_significant_fdr_alpha_0.1_new", mode = "r") # SLP
#         # var_sign= ds_var_sign.variables["masked"][:].data
#         # lons = ds_var.variables["lon"][:].data
#         # lats = ds_var.variables["lat"][:].data
#         lons = np.arange(-180,180,0.5)
#         lats = np.arange(-90,90.5,0.5)

#         # mask_neg = var_sign[:,:] < 0
#         # mask_pos = var_sign[:,:] > 0

#         # ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks/expl_cyclones/bomb_cyclone_freq_{region}_group{k+1}_0_to_3days", mode = "r") # SLP
#         # var_cont = np.array(ds_var.variables["frequency"])

#         var = np.where(var >= 100, np.nan,var)
#         var, lons = add_cyclic_point(var, coord=lons)

#         ax = fig.add_subplot(gs[1,k], projection = map_projection)

#         #create subplots
#         ax.coastlines(color = "black", linewidth = 1.3)
#         ax.set_extent(extent, crs = ccrs.PlateCarree()) 
#         ax.set_facecolor("white")

#         gl = ax.gridlines(draw_labels = True, dms  = True,color = "gray", linewidth = 0.5) # add  gridlines for  latitude and longitude
#         gl.top_labels = False
#         gl.right_labels = False
#         # gl.bottom_labels = False
#         gl.left_labels = False
#         gl.xlabel_style = {"size" : 30, "color" : "k"}
#         gl.ylabel_style = {"size" : 30, "color" : "k"}
#         gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
#         gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
#         # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
#         # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

#         h= ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom ,norm = norm_anom , extend = "both", transform=ccrs.PlateCarree(), alpha =1)
#         i = ax.contour(lons,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
#         # ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
#         # ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())   

#         ax.clabel(i)

        
# cbar_ax = fig.add_axes([0.25, 0.015, 0.5, 0.03])
# cbar = fig.colorbar(h, cax=cbar_ax, orientation= "horizontal")
# cbar.set_label(label = r"frequency anomaly ", size = 30)
# cbar.set_ticks([-0.09,-0.05,-0.01,0.03,0.07])
# cbar.set_ticklabels(["-0.09","-0.05","-0.01","0.03","0.07"])
# cbar.ax.tick_params(labelsize = 30)

# figname = f"/net/thermo/atmosdyn/fschnyde/plots/storm_tracks/all_cyclone_expl_cyclone_freq_anomaly_{region}_3_to_5days.pdf"
# fig.savefig(figname, bbox_inches = "tight",facecolor = "none",edgecolor = "none")
# plt.close(fig)


# %%

# plot_var = "bomb"
# region = "Siberia"
# # region = "NA_sea_ice_allowed"

# titles = ["small BCAR contribution", "mean BCAR contribution", "large BCAR contribution"]

# # cmap_sst
# cmap = matplotlib.cm.get_cmap("binary")
# white = cmap(0)
# cmap = matplotlib.cm.get_cmap("RdBu_r")
# cmap_start = []
# for i in np.linspace(0,0.45,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)
# cmap_start.append(white)
# for i in np.linspace(0.55,1.01,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)    
# cmap_colors = cmap_start
# cmap_sst = ListedColormap(cmap_colors, 'cmap_temp')

# # plot composites
# if region == "Siberia":
#     extent = [120,240,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude=180)
#     xy = (80,40)

# elif region == "NA_sea_ice_allowed":
#     extent = [-120,20,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude= 0)
#     xy = (-90,40)

# cmap_anom = cmap_sst
# levels_anom = np.arange(-0.09,0.11,0.02)
# norm_anom = matplotlib.colors.TwoSlopeNorm(vmin = -0.1, vmax = 0.1, vcenter = 0)
# color_label_anom = r"cyclone frequency anomaly"

# # cmap_sst
# cmap = matplotlib.cm.get_cmap("binary")
# white = cmap(0)
# cmap = matplotlib.cm.get_cmap("RdBu_r")
# cmap_start = []
# for i in np.linspace(0,0.45,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)
# cmap_start.append(white)
# for i in np.linspace(0.55,1.01,100):
#     rgba = cmap(i)
#     cmap_start.append(rgba)    
# cmap_colors = cmap_start
# cmap_sst = ListedColormap(cmap_colors, 'cmap_temp')

# # plot composites
# if region == "Siberia":
#     extent = [120,240,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude=180)
#     xy = (80,40)

# elif region == "NA_sea_ice_allowed":
#     extent = [-100,20,20,90]
#     map_projection = ccrs.PlateCarree(central_longitude= 0)
#     xy = (-90,40)

# cmap_abs = "Oranges"
# levels_abs = np.arange(0,0.5,0.05)
# norm_abs = matplotlib.colors.TwoSlopeNorm(vmin = 0, vmax = 0.5, vcenter = 0.25)
# color_label_abs = r"cyclone frequency anomaly"


# path_clim = "/net/tropo/atmosdyn/mfederer/ape/analysis/clim/cmask_expl_mean_djf_1979_2019"
# clim_ds = nc.Dataset(path_clim, mode = "r")
# clim_expl = np.array(clim_ds.variables["INPUT_SEL"])[0,:,:]

# path_clim = "/net/tropo/atmosdyn/mfederer/ape/analysis/clim/cmask_mean_djf_1979_2019"
# clim_ds = nc.Dataset(path_clim, mode = "r")
# clim_all = np.array(clim_ds.variables["INPUT_SEL"])[0,:,:]

# fig = plt.figure(figsize = (25,18))
# gs = gridspec.GridSpec(3,2, wspace = 0.1)

# # plot weak 
# ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group1_1_to_3days", mode = "r") # SLP
# var = np.array(ds_var.variables["anomaly"])

# if region == "Siberia":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_1_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
# elif region == "NA_sea_ice_allowed":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_NA_group_1_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
# var_sign = ds_var_sign.variables["masked"][:].data

# lons = np.arange(-180,180,0.5)
# lats = np.arange(-90,90.5,0.5)
# lons_clim = np.arange(-180,180,0.5)
# lons_sign = np.arange(-180,180,0.5)
       
# var = np.where(var >= 100, np.nan,var)
# var, lons = add_cyclic_point(var, coord=lons)
# var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

# mask_neg = var_sign[:,:] < 0
# mask_pos = var_sign[:,:] > 0

# clim_expl, lons_clim = add_cyclic_point(clim_expl, coord=lons_clim)

# ax = fig.add_subplot(gs[0,0], projection = map_projection)
# ax.set_title(titles[0], loc = "center", fontsize = 30)

# #create subplots
# ax.coastlines(color = "black", linewidth = 1.3)
# ax.set_extent(extent, crs = ccrs.PlateCarree()) 
# ax.set_facecolor("white")

# gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
# gl.top_labels = False
# gl.right_labels = False
# # gl.bottom_labels = False
# gl.xlabel_style = {"size" : 30, "color" : "k"}
# gl.ylabel_style = {"size" : 30, "color" : "k"}
# gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
# gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))

# h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)           
# i = ax.contour(lons_clim,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
# ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
# ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 

# # plot large BCAR contribution
# ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group3_1_to_3days", mode = "r") # SLP
# var = np.array(ds_var.variables["anomaly"])

# if region == "Siberia":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_3_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
# elif region == "NA_sea_ice_allowed":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_NA_group_3_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
# var_sign= ds_var_sign.variables["masked"][:].data

        
# lons = np.arange(-180,180,0.5)
# lats = np.arange(-90,90.5,0.5)
# lons_sign = np.arange(-180,180,0.5)

# var = np.where(var >= 100, np.nan,var)
# var, lons = add_cyclic_point(var, coord=lons)
# var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

# mask_neg = var_sign[:,:] < 0
# mask_pos = var_sign[:,:] > 0

# ax = fig.add_subplot(gs[0,1], projection = map_projection)
# ax.set_title(titles[2], loc = "center", fontsize = 30)

# #create subplots
# ax.coastlines(color = "black", linewidth = 1.3)
# ax.set_extent(extent, crs = ccrs.PlateCarree()) 
# ax.set_facecolor("white")

# gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
# gl.top_labels = False
# gl.right_labels = False
# gl.left_labels = False
# # gl.bottom_labels = False
# gl.xlabel_style = {"size" : 30, "color" : "k"}
# gl.ylabel_style = {"size" : 30, "color" : "k"}
# gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
# gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

# h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)           
# i = ax.contour(lons_clim,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
# ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
# ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 
# # ax.clabel(i)

# # plot day 3 to 5
# # read in data
# ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group1_3_to_5days", mode = "r") # SLP
# var = np.array(ds_var.variables["anomaly"])

# if region == "Siberia":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_1_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
# elif region == "NA_sea_ice_allowed":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_NA_group_1_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
# var_sign= ds_var_sign.variables["masked"][:].data


# lons = np.arange(-180,180,0.5)
# lats = np.arange(-90,90.5,0.5)
# lons_sign = np.arange(-180,180,0.5)

# var = np.where(var >= 100, np.nan,var)
# var, lons = add_cyclic_point(var, coord=lons)
# var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

# mask_neg = var_sign[:,:] < 0
# mask_pos = var_sign[:,:] > 0

# ax = fig.add_subplot(gs[1,0], projection = map_projection)

# #create subplots
# ax.coastlines(color = "black", linewidth = 1.3)
# ax.set_extent(extent, crs = ccrs.PlateCarree()) 
# ax.set_facecolor("white")

# gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
# gl.top_labels = False
# gl.right_labels = False
# # gl.bottom_labels = False
# # gl.left_labels = False
# gl.xlabel_style = {"size" : 30, "color" : "k"}
# gl.ylabel_style = {"size" : 30, "color" : "k"}
# gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
# gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

# h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)    
# i = ax.contour(lons_clim,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
# ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
# ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 
# # # ax.clabel(i)

# # plot anomaly
# ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group3_3_to_5days", mode = "r") # SLP
# var = np.array(ds_var.variables["anomaly"])

# if region == "Siberia":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_3_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
# elif region == "NA_sea_ice_allowed":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_NA_group_3_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
# var_sign= ds_var_sign.variables["masked"][:].data


# lons = np.arange(-180,180,0.5)
# lats = np.arange(-90,90.5,0.5)
# lons_sign = np.arange(-180,180,0.5)
   
# var = np.where(var >= 100, np.nan,var)
# var, lons = add_cyclic_point(var, coord=lons)
# var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

# mask_neg = var_sign[:,:] < 0
# mask_pos = var_sign[:,:] > 0

# ax = fig.add_subplot(gs[1,1], projection = map_projection)
        
# #create subplots
# ax.coastlines(color = "black", linewidth = 1.3)
# ax.set_extent(extent, crs = ccrs.PlateCarree()) 
# ax.set_facecolor("white")

# gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
# gl.top_labels = False
# gl.right_labels = False
# # gl.bottom_labels = False
# gl.left_labels = False
# gl.xlabel_style = {"size" : 30, "color" : "k"}
# gl.ylabel_style = {"size" : 30, "color" : "k"}
# gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
# gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

# h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom  ,levels = levels_anom ,norm = norm_anom , extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)    
# i = ax.contour(lons_clim,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
# ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
# ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 
# # ax.clabel(i)

# # plot 5 to 7 days
# # read in data
# ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group1_5_to_7days", mode = "r") # SLP
# var = np.array(ds_var.variables["anomaly"])

# if region == "Siberia":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_1_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
# elif region == "NA_sea_ice_allowed":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_NA_group_1_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
# var_sign= ds_var_sign.variables["masked"][:].data

# lons = np.arange(-180,180,0.5)
# lats = np.arange(-90,90.5,0.5)
# lons_sign = np.arange(-180,180,0.5)

# var = np.where(var >= 100, np.nan,var)
# var, lons = add_cyclic_point(var, coord=lons)
# var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)
# mask_neg = var_sign[:,:] < 0
# mask_pos = var_sign[:,:] > 0

# ax = fig.add_subplot(gs[2,0], projection = map_projection)
# # ax.set_title(titles[k], loc = "center", fontsize = 30)

# #create subplots
# ax.coastlines(color = "black", linewidth = 1.3)
# ax.set_extent(extent, crs = ccrs.PlateCarree()) 
# ax.set_facecolor("white")

# gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
# gl.top_labels = False
# gl.right_labels = False
# # gl.bottom_labels = False
# # gl.left_labels = False
# gl.xlabel_style = {"size" : 30, "color" : "k"}
# gl.ylabel_style = {"size" : 30, "color" : "k"}
# gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
# gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

# h= ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)
# i = ax.contour(lons_clim,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())   
# ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
# ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 

# # ax.clabel(i)

# # plot anomaly
# ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group3_5_to_7days", mode = "r") # SLP
# var = np.array(ds_var.variables["anomaly"])

# if region == "Siberia":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_3_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
# elif region == "NA_sea_ice_allowed":
#     ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/expl_cyclones/{plot_var}_cyclone_freq_anom_NA_group_3_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
# var_sign= ds_var_sign.variables["masked"][:].data

# lons = np.arange(-180,180,0.5)
# lats = np.arange(-90,90.5,0.5)
# lons_sign = np.arange(-180,180,0.5)

# var = np.where(var >= 100, np.nan,var)
# var, lons = add_cyclic_point(var, coord=lons)
# var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

# mask_neg = var_sign[:,:] < 0
# mask_pos = var_sign[:,:] > 0

# ax = fig.add_subplot(gs[2,1], projection = map_projection)

# #create subplots
# ax.coastlines(color = "black", linewidth = 1.3)
# ax.set_extent(extent, crs = ccrs.PlateCarree()) 
# ax.set_facecolor("white")

# gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
# gl.top_labels = False
# gl.right_labels = False
# # gl.bottom_labels = False
# gl.left_labels = False
# gl.xlabel_style = {"size" : 30, "color" : "k"}
# gl.ylabel_style = {"size" : 30, "color" : "k"}
# gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
# gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# # gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# # gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

# h= ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom ,norm = norm_anom , extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)
# i = ax.contour(lons_clim,lats, clim_expl[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
# ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
# ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())  

# # ax.clabel(i)

# cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
# cbar = fig.colorbar(h, cax=cbar_ax, orientation= "horizontal")
# cbar.set_label(label = r"frequency anomaly ", size = 30)
# cbar.set_ticks([-0.09,-0.05,-0.01,0.03,0.07])
# cbar.set_ticklabels(["-0.09","-0.05","-0.01","0.03","0.07"])
# cbar.ax.tick_params(labelsize = 30)

# figname = f"/net/thermo/atmosdyn/fschnyde/plots/storm_tracks/{plot_var}_cyclone_freq_anomaly_{region}_1_to_7days.pdf"
# fig.savefig(figname, bbox_inches = "tight",facecolor = "none",edgecolor = "none")
# plt.close(fig)


# %%
plot_var = "all"
region = "Siberia"
# region = "NA_sea_ice_allowed"

titles = ["small BCAR contribution", "mean BCAR contribution", "large BCAR contribution"]

# cmap_sst
cmap = matplotlib.cm.get_cmap("binary")
white = cmap(0)
cmap = matplotlib.cm.get_cmap("RdBu_r")
cmap_start = []
for i in np.linspace(0,0.45,100):
    rgba = cmap(i)
    cmap_start.append(rgba)
cmap_start.append(white)
for i in np.linspace(0.55,1.01,100):
    rgba = cmap(i)
    cmap_start.append(rgba)    
cmap_colors = cmap_start
cmap_sst = ListedColormap(cmap_colors, 'cmap_temp')

# plot composites
if region == "Siberia":
    extent = [120,240,20,90]
    map_projection = ccrs.PlateCarree(central_longitude=180)
    xy = (80,40)

elif region == "NA_sea_ice_allowed":
    extent = [-120,20,20,90]
    map_projection = ccrs.PlateCarree(central_longitude= 0)
    xy = (-90,40)

cmap_anom = cmap_sst
levels_anom = np.arange(-0.09,0.11,0.02)
norm_anom = matplotlib.colors.TwoSlopeNorm(vmin = -0.1, vmax = 0.1, vcenter = 0)
color_label_anom = r"cyclone frequency anomaly"

# cmap_sst
cmap = matplotlib.cm.get_cmap("binary")
white = cmap(0)
cmap = matplotlib.cm.get_cmap("RdBu_r")
cmap_start = []
for i in np.linspace(0,0.45,100):
    rgba = cmap(i)
    cmap_start.append(rgba)
cmap_start.append(white)
for i in np.linspace(0.55,1.01,100):
    rgba = cmap(i)
    cmap_start.append(rgba)    
cmap_colors = cmap_start
cmap_sst = ListedColormap(cmap_colors, 'cmap_temp')

# plot composites
if region == "Siberia":
    extent = [120,240,20,90]
    map_projection = ccrs.PlateCarree(central_longitude=180)
    xy = (80,40)

elif region == "NA_sea_ice_allowed":
    extent = [-100,20,20,90]
    map_projection = ccrs.PlateCarree(central_longitude= 0)
    xy = (-90,40)

cmap_abs = "Oranges"
levels_abs = np.arange(0,0.5,0.05)
norm_abs = matplotlib.colors.TwoSlopeNorm(vmin = 0, vmax = 0.5, vcenter = 0.25)
color_label_abs = r"cyclone frequency anomaly"


path_clim = "/net/tropo/atmosdyn/mfederer/ape/analysis/clim/cmask_expl_mean_djf_1979_2019"
clim_ds = nc.Dataset(path_clim, mode = "r")
clim_expl = np.array(clim_ds.variables["INPUT_SEL"])[0,:,:]

path_clim = "/net/tropo/atmosdyn/mfederer/ape/analysis/clim/cmask_mean_djf_1979_2019"
clim_ds = nc.Dataset(path_clim, mode = "r")
clim_all = np.array(clim_ds.variables["INPUT_SEL"])[0,:,:]

fig = plt.figure(figsize = (25,18))
gs = gridspec.GridSpec(3,2, wspace = 0.1)

# plot weak 
ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group1_1_to_3days", mode = "r") # SLP
var = np.array(ds_var.variables["anomaly"])

if region == "Siberia":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_1_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
elif region == "NA_sea_ice_allowed":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_NA_group_1_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
var_sign = ds_var_sign.variables["masked"][:].data

lons = np.arange(-180,180,0.5)
lats = np.arange(-90,90.5,0.5)
lons_clim = np.arange(-180,180,0.5)
lons_sign = np.arange(-180,180,0.5)
       
var = np.where(var >= 100, np.nan,var)
var, lons = add_cyclic_point(var, coord=lons)
var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

mask_neg = var_sign[:,:] < 0
mask_pos = var_sign[:,:] > 0

clim_all, lons_clim = add_cyclic_point(clim_all, coord=lons_clim)

ax = fig.add_subplot(gs[0,0], projection = map_projection)
ax.set_title(titles[0], loc = "center", fontsize = 30)

#create subplots
ax.coastlines(color = "black", linewidth = 1.3)
ax.set_extent(extent, crs = ccrs.PlateCarree()) 
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
# gl.bottom_labels = False
gl.xlabel_style = {"size" : 30, "color" : "k"}
gl.ylabel_style = {"size" : 30, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))

h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)           
i = ax.contour(lons_clim,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 

# plot large BCAR contribution
ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group3_1_to_3days", mode = "r") # SLP
var = np.array(ds_var.variables["anomaly"])

if region == "Siberia":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_3_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
elif region == "NA_sea_ice_allowed":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_NA_group_3_1_to_3days_significant_fdr_alpha_0.1", mode = "r") # SLP
var_sign= ds_var_sign.variables["masked"][:].data

        
lons = np.arange(-180,180,0.5)
lats = np.arange(-90,90.5,0.5)
lons_sign = np.arange(-180,180,0.5)

var = np.where(var >= 100, np.nan,var)
var, lons = add_cyclic_point(var, coord=lons)
var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

mask_neg = var_sign[:,:] < 0
mask_pos = var_sign[:,:] > 0

ax = fig.add_subplot(gs[0,1], projection = map_projection)
ax.set_title(titles[2], loc = "center", fontsize = 30)

#create subplots
ax.coastlines(color = "black", linewidth = 1.3)
ax.set_extent(extent, crs = ccrs.PlateCarree()) 
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
# gl.bottom_labels = False
gl.xlabel_style = {"size" : 30, "color" : "k"}
gl.ylabel_style = {"size" : 30, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)           
i = ax.contour(lons_clim,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 
# ax.clabel(i)

# plot day 3 to 5
# read in data
ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group1_3_to_5days", mode = "r") # SLP
var = np.array(ds_var.variables["anomaly"])

if region == "Siberia":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_1_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
elif region == "NA_sea_ice_allowed":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_NA_group_1_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
var_sign= ds_var_sign.variables["masked"][:].data


lons = np.arange(-180,180,0.5)
lats = np.arange(-90,90.5,0.5)
lons_sign = np.arange(-180,180,0.5)

var = np.where(var >= 100, np.nan,var)
var, lons = add_cyclic_point(var, coord=lons)
var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

mask_neg = var_sign[:,:] < 0
mask_pos = var_sign[:,:] > 0

ax = fig.add_subplot(gs[1,0], projection = map_projection)

#create subplots
ax.coastlines(color = "black", linewidth = 1.3)
ax.set_extent(extent, crs = ccrs.PlateCarree()) 
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
# gl.bottom_labels = False
# gl.left_labels = False
gl.xlabel_style = {"size" : 30, "color" : "k"}
gl.ylabel_style = {"size" : 30, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)    
i = ax.contour(lons_clim,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 
# # ax.clabel(i)

# plot anomaly
ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group3_3_to_5days", mode = "r") # SLP
var = np.array(ds_var.variables["anomaly"])

if region == "Siberia":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_3_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
elif region == "NA_sea_ice_allowed":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_NA_group_3_3_to_5days_significant_fdr_alpha_0.1", mode = "r") # SLP
var_sign= ds_var_sign.variables["masked"][:].data


lons = np.arange(-180,180,0.5)
lats = np.arange(-90,90.5,0.5)
lons_sign = np.arange(-180,180,0.5)
   
var = np.where(var >= 100, np.nan,var)
var, lons = add_cyclic_point(var, coord=lons)
var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

mask_neg = var_sign[:,:] < 0
mask_pos = var_sign[:,:] > 0

ax = fig.add_subplot(gs[1,1], projection = map_projection)
        
#create subplots
ax.coastlines(color = "black", linewidth = 1.3)
ax.set_extent(extent, crs = ccrs.PlateCarree()) 
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
# gl.bottom_labels = False
gl.left_labels = False
gl.xlabel_style = {"size" : 30, "color" : "k"}
gl.ylabel_style = {"size" : 30, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

h = ax.contourf(lons, lats, var[:,:], cmap = cmap_anom  ,levels = levels_anom ,norm = norm_anom , extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)    
i = ax.contour(lons_clim,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 
# ax.clabel(i)

# plot 5 to 7 days
# read in data
ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group1_5_to_7days", mode = "r") # SLP
var = np.array(ds_var.variables["anomaly"])

if region == "Siberia":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_1_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
elif region == "NA_sea_ice_allowed":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_NA_group_1_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
var_sign= ds_var_sign.variables["masked"][:].data

lons = np.arange(-180,180,0.5)
lats = np.arange(-90,90.5,0.5)
lons_sign = np.arange(-180,180,0.5)

var = np.where(var >= 100, np.nan,var)
var, lons = add_cyclic_point(var, coord=lons)
var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)
mask_neg = var_sign[:,:] < 0
mask_pos = var_sign[:,:] > 0

ax = fig.add_subplot(gs[2,0], projection = map_projection)
# ax.set_title(titles[k], loc = "center", fontsize = 30)

#create subplots
ax.coastlines(color = "black", linewidth = 1.3)
ax.set_extent(extent, crs = ccrs.PlateCarree()) 
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
# gl.bottom_labels = False
# gl.left_labels = False
gl.xlabel_style = {"size" : 30, "color" : "k"}
gl.ylabel_style = {"size" : 30, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

h= ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom,norm = norm_anom, extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)
i = ax.contour(lons_clim,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())   
ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree()) 

# ax.clabel(i)

# plot anomaly
ds_var = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anomaly_{region}_group3_5_to_7days", mode = "r") # SLP
var = np.array(ds_var.variables["anomaly"])

if region == "Siberia":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_Siberia_group_3_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
elif region == "NA_sea_ice_allowed":
    ds_var_sign = nc.Dataset(f"/net/thermo/atmosdyn/fschnyde/data/storm_tracks_final/cyclones/{plot_var}_cyclone_freq_anom_NA_group_3_5_to_7days_significant_fdr_alpha_0.1", mode = "r") # SLP
var_sign= ds_var_sign.variables["masked"][:].data

lons = np.arange(-180,180,0.5)
lats = np.arange(-90,90.5,0.5)
lons_sign = np.arange(-180,180,0.5)

var = np.where(var >= 100, np.nan,var)
var, lons = add_cyclic_point(var, coord=lons)
var_sign, lons_sign = add_cyclic_point(var_sign, coord=lons_sign)

mask_neg = var_sign[:,:] < 0
mask_pos = var_sign[:,:] > 0

ax = fig.add_subplot(gs[2,1], projection = map_projection)

#create subplots
ax.coastlines(color = "black", linewidth = 1.3)
ax.set_extent(extent, crs = ccrs.PlateCarree()) 
ax.set_facecolor("white")

gl = ax.gridlines(draw_labels = True, dms  = True,color = "none", linewidth = 0.5) # add  gridlines for  latitude and longitude
gl.top_labels = False
gl.right_labels = False
# gl.bottom_labels = False
gl.left_labels = False
gl.xlabel_style = {"size" : 30, "color" : "k"}
gl.ylabel_style = {"size" : 30, "color" : "k"}
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))
gl.ylocator = mticker.FixedLocator(np.arange(0,100,20))
# gl.xlocator = mticker.FixedLocator([80,100,120,140,160,-180,-160,-140,-120,-100,-80])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70,80])

h= ax.contourf(lons, lats, var[:,:], cmap = cmap_anom ,levels = levels_anom ,norm = norm_anom , extend = "both", transform=ccrs.PlateCarree(), alpha =0.9)
i = ax.contour(lons_clim,lats, clim_all[:,:], levels = np.arange(0.1,0.45,0.05), colors = "gray", linewidths = 2,transform = ccrs.PlateCarree())
ax.contourf(lons, lats, mask_neg, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())     
ax.contourf(lons, lats, mask_pos, hatches=['.'],colors = "none", levels = [0.5,1],transform = ccrs.PlateCarree())  

# ax.clabel(i)

cbar_ax = fig.add_axes([0.2, 0.05, 0.6, 0.02])
cbar = fig.colorbar(h, cax=cbar_ax, orientation= "horizontal")
cbar.set_label(label = r"frequency anomaly ", size = 30)
cbar.set_ticks([-0.09,-0.05,-0.01,0.03,0.07])
cbar.set_ticklabels(["-0.09","-0.05","-0.01","0.03","0.07"])
cbar.ax.tick_params(labelsize = 30)

figname = f"/net/thermo/atmosdyn/fschnyde/plots/storm_tracks/{plot_var}_cyclone_freq_anomaly_{region}_1_to_7days.png"
fig.savefig(figname, bbox_inches = "tight",facecolor = "none",edgecolor = "none", dpi = 300)
plt.close(fig)
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
import scipy.stats
from matplotlib.patches import Rectangle
import operator
from dypy.small_tools  import interpolate

def _get_segments(self, trajs):
    lat, lon = trajs["lat"].squeeze(), trajs["lon"].squeeze()
    points = np.array([lon, lat]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return segments


def plot_traj_varsel(self, data, ax, ind_sel, variable="p"):
    # dimx_lon are the trajectory indices

    segments, colors = list(
        zip(
            *[
                (
                    self._get_segments(
                        data.tra.sel(dimx_lon=i, time=slice(ind_sel[0], ind_sel[1]))
                    ),
                    data.tra.sel(
                        dimx_lon=i,
                        time=slice(
                            ind_sel[0], ind_sel[1] - 1
                        ),  # colors needs to be shorter by 1 time step compared to segments
                    )[variable].values.squeeze(),
                )
                for i in data.itra_sel  # traj indices
            ]
        )
    )

    segments = np.concatenate(segments)
    colors = np.concatenate(colors)
    cmap = matplotlib.cm.get_cmap("RdBu").reversed()

    if variable == "p":
        tmp = data.tra[variable].values.squeeze()
        minlev = np.nanmin(tmp[tmp > 0])
        maxlev = np.nanmax(tmp[tmp > 0])
        levels = np.linspace(minlev, maxlev, 21)

    norm = matplotlib.colors.BoundaryNorm(levels, cmap.N)

    nkwargs = {
        "array": colors,
        "cmap": cmap,
        "norm": norm,
        "linewidths": 0.5,
        "transform": ccrs.PlateCarree(),
    }

    lc = matplotlib.collections.LineCollection(segments, **nkwargs)

    ax.add_collection(lc)

    return lc


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

def longest(numbers):
    max, count_ = 1, 1
    start_idx, end_idx = 0, 0
    for i in range(len(numbers)-1):
        # if difference between number and his follower is 1,they are in sequence
        if numbers[i+1]-numbers[i] ==1:
            count_ = count_+1
        else:
            if count_ > max :
                max = count_
                end_idx = i
                start_idx = i+1 - max
            # Reset counter
            count_ = 1
    return (start_idx,end_idx,max)

plt.style.use('default') # I personally prefer seaborn for the graph style, but you may choose whichever you want.
params = {"ytick.color" : "black",
          "xtick.color" : "black",
          "grid.color" :  (217/256, 217/256, 217/256), #(217/256, 217/256, 217/256), "slategrey"
          "axes.labelcolor" : "black",
          "axes.edgecolor" : "black",
          "text.usetex" : True,
          "font.family" : "serif",
          "font.weight" :"bold",
          "font.serif" : ["Computer Modern Serif"]}
plt.rcParams.update(params)
# %%
# # ----------------------
# # find areas of cooling
# # -----------------------
# integrated over entire trajectory (each point weighted by TH change)
path = "/net/thermo/atmosdyn/fschnyde/data/lagranto/clim_final/Siberia"
tras_BCAR = Tra(path + "/CAO_BCAR_trajs_1979:2022")

lons = np.zeros(len(tras_BCAR)*(len(tras_BCAR[0])-2))
lats = np.zeros(len(tras_BCAR)*(len(tras_BCAR[0])-2))
weights = np.zeros(len(tras_BCAR)*(len(tras_BCAR[0])-2))
counter = 0
for i in range(len(tras_BCAR)):
    tra = tras_BCAR[i].data[::-1]           # get trajectory
    for j in range(6,len(tras_BCAR[0])-6):
        tra_t = tra[j]                      # get lon and lat at specific time of trajectory (for each timestep)
        delta_theta = (tra["TH"][j-6] - tra["TH"][j+6])/12 # get change of Theta in location of trajectory (3h rate), this will be our weight
        if delta_theta > 0:
            weights[counter] = delta_theta
            lon_t = tra_t["lon"]
            lat_t = tra_t["lat"]
            lons[counter] = lon_t
            lats[counter] = lat_t
            counter += 1
        else:
            continue
        
# # %%
df = pd.DataFrame(list(zip(lons, lats, weights)), columns = ["lons","lats","weights"])
# df_new = df.loc[~(df==0).all(axis=1)]
# filename = "/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate"
# df_new.to_csv(filename, sep = " ", index = False, header = False) 

# # %%
lons = np.zeros(len(tras_BCAR)*(len(tras_BCAR[0])-2))
lats = np.zeros(len(tras_BCAR)*(len(tras_BCAR[0])-2))
weights = np.zeros(len(tras_BCAR)*(len(tras_BCAR[0])-2))
counter = 0
for i in range(len(tras_BCAR)):
    tra = tras_BCAR[i].data[::-1]           # get trajectory
    for j in range(6,len(tras_BCAR[0])-6):
        tra_t = tra[j]                      # get lon and lat at specific time of trajectory (for each timestep)
        delta_theta = (tra["TH"][j-6] - tra["TH"][j+6])/12 # get change of Theta in location of trajectory (3h rate), this will be our weight
        if delta_theta > 0:
            weights[counter] = 1
            lon_t = tra_t["lon"]
            lat_t = tra_t["lat"]
            lons[counter] = lon_t
            lats[counter] = lat_t
            counter += 1
        else:
            continue
# # %%
# df = pd.DataFrame(list(zip(lons, lats, weights)), columns = ["lons","lats","weights"])
# df_new = df.loc[~(df==0).all(axis=1)]
# filename = "/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate"
# df_new.to_csv(filename, sep = " ", index = False, header = False) 

# %%
# ----------------------
# find areas of cooling based on 3 categories
# -----------------------
# integrated over entire trajectory (each point weighted by TH change)

for region in ["Siberia"]:
    if region == "Siberia":
        event_path = "/net/thermo/atmosdyn/fschnyde/data/cao/"
        events = pd.read_csv(event_path + "dates_CAO_max_intens_1979_2022_Siberia_final", sep = " ",header = None,
                            names = ["date", "intensity", "rank", "unknown"])

        path = "/net/thermo/atmosdyn/fschnyde/data/lagranto/clim_final/Siberia"
        cao_contributions = pd.read_csv(path +  "/CAO_contributions_1979:2022", sep = ",", skiprows = 1,
                                            names = ["date", "BCAR", "DI", "ul", "ll", "TH_ocean"])
        wcb_contributions = pd.read_csv(path +  "/WCB_contributions_1979:2022", sep = ",", skiprows = 1,
                                            names = ["date", "WCB_3d", "WCB_anytime", "non_BCAR_WCB", "ascent_time","TH_ocean", "microwaving_WCB", "microwaving_nonWCB" ])

        merge_merge = cao_contributions.merge(wcb_contributions, left_on="date", right_on="date")
        merge_merge.TH_ocean_y.fillna(merge_merge.TH_ocean_x, inplace = True)

        q1 = np.nanquantile(merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())), 0.33)
        q2 = np.nanquantile(merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())), 0.66)
        group1 = merge_merge[merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) < q1]
        group2 = merge_merge[(merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) > q1)
                        & (merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) < q2)]
        group3 = merge_merge[merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) > q2]

        groups = [group1,group2,group3]
    elif region == "NA_sea_ice_allowed":
        event_path = "/net/thermo/atmosdyn/fschnyde/data/cao/"
        events = pd.read_csv(event_path + "dates_CAO_max_intens_1979_2022_NA_final", sep = " ",header = None,
                            names = ["date", "intensity", "rank", "unknown"])

        path = "/net/thermo/atmosdyn/fschnyde/data/lagranto/clim_final/NA"
        cao_contributions = pd.read_csv(path +  "/CAO_contributions_1979:2022_NA_sea_ice_allowed", sep = ",", skiprows = 1,
                                            names = ["date", "BCAR", "DI", "ul", "ll", "TH_ocean"])
        wcb_contributions = pd.read_csv(path +  "/WCB_contributions_1979:2022_NA", sep = ",", skiprows = 1,
                                            names = ["date", "WCB_3d", "WCB_anytime", "non_BCAR_WCB", "ascent_time","TH_ocean", "microwaving_WCB", "microwaving_nonWCB" ])

        merge_merge = cao_contributions.merge(wcb_contributions, left_on="date", right_on="date")
        merge_merge.TH_ocean_y.fillna(merge_merge.TH_ocean_x, inplace = True)

        q1 = np.nanquantile(merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())), 0.33)
        q2 = np.nanquantile(merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())), 0.66)
        group1 = merge_merge[merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) < q1]
        group2 = merge_merge[(merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) > q1)
                        & (merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) < q2)]
        group3 = merge_merge[merge_merge["BCAR"]*(1-(merge_merge["TH_ocean_y"]-merge_merge["TH_ocean_y"].min())/(merge_merge["TH_ocean_y"].max()-merge_merge["TH_ocean_y"].min())) > q2]

        groups = [group1,group2,group3]

    # for i in range(3):
    i = 2
    group = groups[i]
    print(group)

    lons = np.zeros(330000*(167))
    lats = np.zeros(330000*(167))
    weights = np.zeros(330000*(167))
    counter = 0

    lons_2days = np.zeros(330000*(167))
    lats_2days = np.zeros(330000*(167))
    weights_2days = np.zeros(330000*(167))
    counter_2days = 0

    lons_4days = np.zeros(330000*(167))
    lats_4days = np.zeros(330000*(167))
    weights_4days = np.zeros(330000*(167))
    counter_4days = 0

    lons_6days = np.zeros(330000*(167))
    lats_6days = np.zeros(330000*(167))
    weights_6days = np.zeros(330000*(167))
    counter_6days = 0


    for j in range(len(group)):

        # go through each CAO event
        date = group["date"].iloc[j]
        print(date)
        if region == "Siberia": 
            path = "/net/thermo/atmosdyn/fschnyde/data/lagranto/clim_final/Siberia/CAO_BCAR_trajs"
            tras_BCAR = Tra(path + f"/CAO_BCAR_trajs_{date}")
        elif region == "NA_sea_ice_allowed": 
            path = "/net/thermo/atmosdyn/fschnyde/data/lagranto/clim_final/NA/CAO_BCAR_trajs"
            tras_BCAR = Tra(path + f"/CAO_BCAR_trajs_{date}_NA_sea_ice_allowed")

        
        for k in range(len(tras_BCAR)):
            tra = tras_BCAR[k].data[::-1]           # get trajectory
            for l in range(6,len(tras_BCAR[0])-6):
                tra_t = tra[l]                      # get lon and lat at specific time of trajectory (for each timestep)
                delta_theta = (tra["TH"][l-6] - tra["TH"][l+6])/12 # get change of Theta in location of trajectory (3h rate), this will be our weight
                if delta_theta > 0:
                    weights[counter] = delta_theta
                    lon_t = tra_t["lon"]
                    lat_t = tra_t["lat"]
                    lons[counter] = lon_t
                    lats[counter] = lat_t
                    counter += 1
                else:
                    continue

                # additionally add to 2days/4days/6days lists

                if l > 108 and l < 132:
                    if delta_theta > 0:
                        weights_2days[counter_2days] = delta_theta
                        lon_t = tra_t["lon"]
                        lat_t = tra_t["lat"]
                        lons_2days[counter_2days] = lon_t
                        lats_2days[counter_2days] = lat_t
                        counter_2days += 1
                    else:
                        continue

                if l > 60 and l < 84:
                    if delta_theta > 0:
                        weights_4days[counter_4days] = delta_theta
                        lon_t = tra_t["lon"]
                        lat_t = tra_t["lat"]
                        lons_4days[counter_4days] = lon_t
                        lats_4days[counter_4days] = lat_t
                        counter_4days += 1
                    else:
                        continue

                if l > 12 and l < 36:
                    if delta_theta > 0:
                        weights_6days[counter_6days] = delta_theta
                        lon_t = tra_t["lon"]
                        lat_t = tra_t["lat"]
                        lons_6days[counter_6days] = lon_t
                        lats_6days[counter_6days] = lat_t
                        counter_6days += 1
                    else:
                        continue

        weights_ones = np.ones(len(weights))

        print(str(group) + " done")

        df = pd.DataFrame(list(zip(lons, lats, weights)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate_group{i}"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_NA_12h_cooling_rate_group{i}"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 

        df = pd.DataFrame(list(zip(lons, lats, weights_ones)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group{i}"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group{i}"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 

        # save 2 day list

        weights_ones_2days = np.ones(len(weights_2days))
        
        df = pd.DataFrame(list(zip(lons_2days, lats_2days, weights_2days)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate_group{i}_2days"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_NA_12h_cooling_rate_group{i}_2days"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 


        df = pd.DataFrame(list(zip(lons_2days, lats_2days, weights_ones_2days)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group{i}_2days"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group{i}_2days"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 

        # save 4 day list

        weights_ones_4days = np.ones(len(weights_4days))
        
        df = pd.DataFrame(list(zip(lons_4days, lats_4days, weights_4days)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate_group{i}_4days"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_NA_12h_cooling_rate_group{i}_4days"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 


        df = pd.DataFrame(list(zip(lons_4days, lats_4days, weights_ones_4days)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group{i}_4days"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group{i}_4days"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 

        # save 6 day list

        weights_ones_6days = np.ones(len(weights_6days))
        
        df = pd.DataFrame(list(zip(lons_6days, lats_6days, weights_6days)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate_group{i}_6days"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_cooling_along_trajectory_DTHweighted_NA_12h_cooling_rate_group{i}_6days"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 


        df = pd.DataFrame(list(zip(lons_6days, lats_6days, weights_ones_6days)), columns = ["lons","lats","weights"])
        df_new = df.loc[~(df==0).any(axis=1)]
        if region == "Siberia":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group{i}_6days"
        elif region == "NA_sea_ice_allowed":
            filename = f"/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group{i}_6days"
        df_new.to_csv(filename, sep = " ", index = False, header = False) 



# %%
df = pd.read_csv("/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_6days", sep = " ",header = None,
                            names = ["lon", "lat", "weight"])
df_new = df.loc[~(df==0).any(axis=1)]
df_new.to_csv("/net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final/gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_6days", sep = " ", index = False, header = False) 
# %%

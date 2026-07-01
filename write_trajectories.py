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
import glob
import os
import itertools



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
member_ids = [1000, 1100, 1200, 1300, 1400]
regions = ['JS', 'NA']
timeframes = ['present', 'future']

# member_ids = [1000, 1100, 1200, 1300, 1400]
# regions = ['NA']
# timeframes = ['present']


# %%
# Identify different air streams in trajectories
# CAO_ALL_TRAJS_..._.csv combines the 4 different airstreams -> used for statistics
# CAO_ALL_TRAJS_..._.txt has all trajectories in ASCII format -> used for plotting
for region in regions:
    for timeframe in timeframes:
        all_all = []
        BCAR_all = []
        DI_all = []
        non_BCAR_all = []
        ll_all = []
        ul_all = []

        # BCAR_all = Tra()
        # DI_all = Tra()
        # non_BCAR_all = Tra()
        # ll_all = Tra()

        if timeframe == 'present':
            Bfile = "B_19*_with_TH.txt"
            Bfile2 = "B_2000*_with_TH.txt"
        elif timeframe == 'future':
            Bfile = "B_209*_with_TH.txt"
            Bfile2 = "B_21*_with_TH.txt"

        for member_id in member_ids:
            path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/lagranto/traj/{region}/"
            files = glob.glob(os.path.join(path, Bfile)) + glob.glob(os.path.join(path, Bfile2))

            for file in files:
                print(f"Processing {file}")
                all_traj = pd.read_csv(file)                
                try:
                    trajs = Tra(file)
                except OSError as e:
                    print(f"Skipping {file}: {e}")
                    continue
                # --- Add this block to check for TH_1 ---
                if 'TH_1' in trajs.variables:
                    print(f"File {file} contains variable TH_1! Variables: {trajs.variables}")
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
                if region == 'JS':
                    cold_inds = np.where((low_trajs["TH"]<280) & (low_trajs["SST"]==0))
                elif region == 'NA':
                    # ll colder than 280K and over land
                    cold_inds = np.where(((low_trajs["TH"]<280)
                    & (low_trajs["SST"]==0))
                    | 

                    # trajectories colder than 280K and over sea-ice in certain box
                    (low_trajs["TH"]<280)
                    & (low_trajs["SST"]<273)
                    & 

                    (low_trajs["lon"]>240)
                    & (low_trajs["lon"]<295)
                    &

                    (low_trajs["lat"]>40)
                    & (low_trajs["lat"]<75))
                
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

                # select that at least at one timestep p has to be higher than 1020hPa
                anti_inds = np.where(cold_trajs["PSL"]>1020)
                anti_inds = np.unique(anti_inds[0])
                anti_trajs = Tra()
                anti_trajs.set_array(cold_trajs[anti_inds,:])

                BCAR_trajs = anti_trajs

                # # OLD: nonBCAR = low level - cold trajectories
                # #  remove non cold trajs from array for non BCAR trajetories
                # if len(cold_trajs) != len(low_trajs):
                #     non_cold = np.delete(array, inds)
                # else: 
                #     non_cold = array

                # if len(cold_trajs) != len(anti_trajs):
                #     print("no anti inds " + file)

                # non_BCAR_trajs = Tra()
                # non_BCAR_trajs.set_array(low_trajs[non_cold,:])

                # NEW: nonBCAR = low level - BCAR trajectories
                bcar_indices_in_low = [inds[i] for i in anti_inds]  # map BCAR indices back to low_trajs
                if len(bcar_indices_in_low) != 0:
                    non_bcar_indices = np.delete(array, bcar_indices_in_low)  # removes BCAR trajectories
                else:
                    non_bcar_indices = array
                    print("no bcar inds " + file)
                non_BCAR_trajs = Tra()
                non_BCAR_trajs.set_array(low_trajs[non_bcar_indices,:])



                # fraction of low level trajectories
                ll_inds = np.where(trajs["p"]>700)
                inds = []
                for j in range(len(trajs)):
                    count  = len(np.where(ll_inds[0] == j)[0])     # get those indices where condition is met for all timesteps
                    if count == 29:
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
                    if count == 29:
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

                if len(BCAR_trajs) > 0:
                    BCAR_all.append(BCAR_trajs)
                if len(DI_tras_good) > 0:
                    DI_all.append(DI_tras_good)
                if len(non_BCAR_trajs) > 0:
                    non_BCAR_all.append(non_BCAR_trajs)
                if len(ll_trajs) > 0:
                    ll_all.append(ll_trajs)
                if len(ul_trajs) > 0:
                    ul_all.append(ul_trajs)
                if len(trajs) > 0:
                    all_all.append(trajs)

                print(f"Number of LL trajectories: {len(ll_trajs)}")
                print(f"Number of BCAR trajectories: {len(BCAR_trajs)}")
                print(f'Number of non-BCAR trajectories: {len(non_BCAR_trajs)}')
                if len(ll_trajs) == len(BCAR_trajs):
                # Check if all trajectories are identical
                    identical = np.all(ll_trajs[:] == BCAR_trajs[:])
                    print(f"LL and BCAR trajectories are identical: {identical}")
                else:
                    print("LL and BCAR have different numbers of trajectories.")
                print(f"Number of UL trajectories: {len(ul_trajs)}")

                # if len(BCAR_trajs) > 0:
                #     BCAR_all = BCAR_all.concatenate(BCAR_trajs)
                # if len(DI_tras_good) > 0:
                #     DI_all = DI_all.concatenate(DI_tras_good)
                # if len(non_BCAR_trajs) > 0:
                #     non_BCAR_all = non_BCAR_all.concatenate(non_BCAR_trajs)
                # if len(ll_trajs) > 0:
                #     ll_all = ll_all.concatenate(ll_trajs)
        outdir = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/"
       
        # Concatenate all trajectories for this region and write to ASCII
        def concat_and_write(traj_list, outname):
            if traj_list:
                first_traj = traj_list[0]
                other_trajs = traj_list[1:]
                # for t in traj_list[1:]:
                #     all_trajs.concatenate(t)
                # all_trajs.write_ascii(outname)
                for i, t in enumerate(traj_list):
                    print(f"Traj {i} variables: {t.variables}")
                all_trajs = first_traj.concatenate(other_trajs)
                all_trajs.write_ascii(outname)

        concat_and_write(BCAR_all, f"{outdir}CAO_BCAR_trajs_{region}_{timeframe}_ALL.txt")
        concat_and_write(DI_all, f"{outdir}CAO_DI_trajs_{region}_{timeframe}_ALL.txt")
        concat_and_write(non_BCAR_all, f"{outdir}CAO_nonBCAR_trajs_{region}_{timeframe}_ALL.txt")
        concat_and_write(ll_all, f"{outdir}CAO_LL_trajs_{region}_{timeframe}_ALL.txt")
        concat_and_write(ul_all, f"{outdir}CAO_UL_trajs_{region}_{timeframe}_ALL.txt")
        concat_and_write(all_all, f"{outdir}CAO_ALL_trajs_{region}_{timeframe}_ALL.txt")
        
        
        def tra_to_df(tra, airstream, member_id, region, traj_offset=0):
            """Convert a Tra object to a DataFrame with metadata columns."""
            if len(tra) == 0:
                return None
            arr = tra[:]  # This is a structured array: shape (n_traj, n_step)
            n_traj, n_step = arr.shape[0], arr.shape[1]
            arr_flat = arr.reshape(-1)  # shape (n_traj * n_step,)
            df = pd.DataFrame({name: arr_flat[name] for name in arr.dtype.names})
            df['traj_id'] = np.repeat(np.arange(n_traj) + traj_offset, n_step)
            df['airstream'] = airstream
            df['member_id'] = member_id

            # --- Convert time to relative hours per trajectory ---
            # If 'time' is datetime, convert to relative hours
            if np.issubdtype(df['time'].dtype, np.datetime64):
                df['rel_time'] = (
                    df.groupby('traj_id')['time']
                    .transform(lambda x: (x - x.iloc[0]).dt.total_seconds() / 3600)
                    .round(2)
                )
            # If 'time' is string, try to convert to datetime first
            elif df['time'].dtype == object:
                try:
                    df['time'] = pd.to_datetime(df['time'])
                    df['rel_time'] = (
                        df.groupby('traj_id')['time']
                        .transform(lambda x: (x - x.iloc[0]).dt.total_seconds() / 3600)
                        .round(2)
                    )
                except Exception:
                    # If conversion fails, leave as is
                    df['rel_time'] = df['time']
            else:
                # If already numeric, set relative to first value in each trajectory
                df['rel_time'] = (
                    df.groupby('traj_id')['time']
                    .transform(lambda x: x - x.iloc[0])
                    .round(2)
                )
            return df

        # Collect all DataFrames for this region/timeframe
        dfs = []
        traj_offset = 0

        for tra_list, airstream in zip([BCAR_all, DI_all, non_BCAR_all, ll_all, ul_all],
                                    ['BCAR', 'DI', 'nonBCAR', 'LL', 'UL']):
            for i, tra in enumerate(tra_list):
                df = tra_to_df(tra, airstream, member_ids[0], region, traj_offset)
                if df is not None:
                    dfs.append(df)
                    traj_offset += len(tra)

        if dfs:
            all_trajs_df = pd.concat(dfs, ignore_index=True)
            # Save as CSV for easy plotting
            out_csv = f"{outdir}CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
            all_trajs_df.to_csv(out_csv, index=False)
            print(f"Saved combined CSV: {out_csv}")
        # BCAR_all.write_ascii(f"{outdir}CAO_BCAR_trajs_{region}_{timeframe}_ALL.txt")
        # DI_all.write_ascii(f"{outdir}CAO_DI_trajs_{region}_{timeframe}_ALL.txt")
        # non_BCAR_all.write_ascii(f"{outdir}CAO_nonBCAR_low_level_trajs_{region}_{timeframe}_ALL.txt")
        # ll_all.write_ascii(f"{outdir}CAO_all_low_level_trajs_{region}_{timeframe}_ALL.txt")
# %%
# txt version with all trajectories combined for plotting trajectories
for region in regions:
    for timeframe in timeframes:
        all_all = []
        BCAR_all = []
        DI_all = []
        non_BCAR_all = []
        ll_all = []

        # BCAR_all = Tra()
        # DI_all = Tra()
        # non_BCAR_all = Tra()
        # ll_all = Tra()

        if timeframe == 'present':
            Bfile = "B_19*_with_TH.txt"
            Bfile2 = "B_2000*_with_TH.txt"
        elif timeframe == 'future':
            Bfile = "B_209*_with_TH.txt"
            Bfile2 = "B_21*_with_TH.txt"

        for member_id in member_ids:
            path = f"/net/litho/atmosdyn2/nfieldhouse/data/{member_id}/lagranto/traj/{region}/"
            files = glob.glob(os.path.join(path, Bfile)) + glob.glob(os.path.join(path, Bfile2))

            for file in files:
                print(f"Processing {file}")
                all_traj = pd.read_csv(file)                
                try:
                    trajs = Tra(file)
                except OSError as e:
                    print(f"Skipping {file}: {e}")
                    continue
                # --- Add this block to check for TH_1 ---
                if 'TH_1' in trajs.variables:
                    print(f"File {file} contains variable TH_1! Variables: {trajs.variables}")
                            # -------------------------------
                #  identify different air streams
                # -------------------------------
                ###### select BCAR trajectories
                # select that p always larger than 700hPa

                all_all.append(trajs)

                # if len(BCAR_trajs) > 0:
                #     BCAR_all = BCAR_all.concatenate(BCAR_trajs)
                # if len(DI_tras_good) > 0:
                #     DI_all = DI_all.concatenate(DI_tras_good)
                # if len(non_BCAR_trajs) > 0:
                #     non_BCAR_all = non_BCAR_all.concatenate(non_BCAR_trajs)
                # if len(ll_trajs) > 0:
                #     ll_all = ll_all.concatenate(ll_trajs)
        outdir = "/net/litho/atmosdyn2/nfieldhouse/data/lagranto/"
       
        # Concatenate all trajectories for this region and write to ASCII
        def concat_and_write(traj_list, outname):
            if traj_list:
                first_traj = traj_list[0]
                other_trajs = traj_list[1:]
                # for t in traj_list[1:]:
                #     all_trajs.concatenate(t)
                # all_trajs.write_ascii(outname)
                for i, t in enumerate(traj_list):
                    print(f"Traj {i} variables: {t.variables}")
                all_trajs = first_traj.concatenate(other_trajs)
                all_trajs.write_ascii(outname)
        
        concat_and_write(all_all, f"{outdir}CAO_ALL_trajs_{region}_{timeframe}_ALL.txt")
# %%

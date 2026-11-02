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
import xarray as xr
import scipy.stats
import matplotlib.dates as mdates  # Import mdates for date formatting
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
import config

# %%
member_id=sys.argv[1]  # Get from command line
base_path = f"{config.DATA_ROOT}/{member_id}/cao/clim"
fig_path = f"{config.PLOTS_ROOT}/"
regions = ["JS", "NA"]
timeframes = ["present", "future"]

for region in regions:
    for timeframe in timeframes:
        path = f"{base_path}/{region}/"
        ds_anom = xr.open_dataset(path + f"yearly_hour_anomalies_{timeframe}", decode_times=True, engine='netcdf4')

        # old method
        # CAOI = ds_anom.variables["CAOI"].values
        # CAOI_mean = np.nanmean(np.nanmean(CAOI, axis=2), axis=2)[:, 0]
        # # CAOI_min = np.nanmin(np.nanmin(CAOI, axis=2), axis=2)[:, 0]

        # use latitudinal weights to calculate the mean
        CAOI = ds_anom.CAOI
        weights = np.cos(np.deg2rad(CAOI.lat))
        weights.name = "weights"
        CAOI_weighted = CAOI.weighted(weights)
        CAOI_mean = CAOI_weighted.mean(dim=["lat", "lon", "lev"], skipna=True).values

        # make list (will be input for the events script)
        dates = ds_anom.indexes["time"].to_datetimeindex()  # set index time of the dataarray to datetime
        dates_fmt = []
        for date in dates:
            date_str = dt.datetime.strftime(date, '%Y-%m-%d-%H')
            dates_fmt.append(date_str)
        data_dict = {"dates": dates_fmt, "CAOI": CAOI_mean}
        df = pd.DataFrame(data_dict)

        fout = f"{path}CAOI_{timeframe}.csv"
        df.to_csv(fout, sep=" ", header=None, index=False)

        # Select only DJF months
        ds_anom_djf = ds_anom.sel(time=ds_anom['time.season'] == 'DJF')
        
        # read variables
        CAOI_djf = ds_anom_djf.CAOI
        CAOI_djf_weighted = CAOI_djf.weighted(weights)
        CAOI_djf_mean = CAOI_djf_weighted.mean(dim=["lat", "lon", "lev"], skipna=True)
        
        # plot time series of CAOI mean (and plot threshold for 10th percentile --> strongest CAOs)
        q10 = np.quantile(CAOI_djf_mean, 0.1)
        plt.figure(figsize=(10, 5))
        plt.plot(ds_anom.indexes["time"].to_datetimeindex(), CAOI_mean, label='Mean')  # Use dates for x-axis
        plt.axhline(q10, color="red", linestyle='--', label='10th Percentile of DJF')
        plt.title(f"CAOI Anomaly Mean within Box for {region} {timeframe}")
        plt.xlabel('Time')
        plt.ylabel('CAOI Anomaly')
        plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
        plt.legend(loc='upper right')
        plt.grid(True)

        # Format the x-axis to show fewer dates
        ax = plt.gca()
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))  # Show every 6 months
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))  # Format as Year-Month

        # Adjust layout to prevent cutting off labels
        plt.tight_layout()

        # Save the figure
        fig_filename = f"{fig_path}{member_id}_timeseries_CAOIano_{region}_{timeframe}.png"
        plt.savefig(fig_filename)
        # plt.show()
        plt.close()  # Close the figure to free up memory


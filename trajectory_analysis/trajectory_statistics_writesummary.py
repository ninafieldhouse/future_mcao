# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import sys
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
import config

plt.rcParams.update({
    'font.size': 16,           # Default font size for ticks
    'axes.titlesize': 20,      # Title font size
    'axes.labelsize': 20,      # X and Y label font size
    'legend.fontsize': 20,     # Legend font size
    'xtick.labelsize': 16,     # X-axis tick font size
    'ytick.labelsize': 16      # Y-axis tick font size
})

# plt.rcParams.update({
#     'font.size': 8,           # Default font size for ticks
#     'axes.titlesize': 8,      # Title font size
#     'axes.labelsize': 8,      # X and Y label font size
#     'legend.fontsize': 6,     # Legend font size
#     'xtick.labelsize': 6,     # X-axis tick font size
#     'ytick.labelsize': 6      # Y-axis tick font size
# })

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"

pt = 1./72.27 # Hundreds of years of history... 72.27 points to an inch.

textwidth = 418.25368


figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
            # CQG is only one column
            # Add more journals below. Can add more properties to each journal

my_width = figsizes["twocol"]
# Our figure's aspect ratio
golden = (1 + 5 ** 0.5) / 2
fourto3 = 4/3

regions = ['JS', 'NA']
airstreams = ['LL', 'UL']
timeframes = ['present', 'future']
colors = {'present': 'blue', 'future': 'orange'}
figpath = f'{config.PLOTS_ROOT}/traj/'
datapath = f'{config.LAGRANTO_DIR}/'

units_dict = {
    'Q': 'g kg${-1}$',
    'TH': 'K',
    'SHFLX': 'W m${-2}$',
    'LHFLX': 'W m${-2}$'
}

def weighted_mean(values, weights):
    return np.sum(values * weights) / np.sum(weights)

def weighted_quantile(values, weights, quantile):
    sorter = np.argsort(values)
    values, weights = values[sorter], weights[sorter]
    cumulative_weight = np.cumsum(weights)
    cutoff = quantile * cumulative_weight[-1]
    return values[np.searchsorted(cumulative_weight, cutoff)]

# def calculate_weighted_stats(df, varname):
#     """Returns dict: rel_time -> (mean, q25, q75) for the given variable, latitude-weighted."""
#     df = df.copy()
#     df['weight'] = np.cos(np.deg2rad(df['lat']))
#     stats = {}
#     for t, group in df.groupby('rel_time'):
#         vals = group[varname].values
#         w = group['weight'].values
#         stats[t] = (
#             weighted_mean(vals, w),
#             weighted_quantile(vals, w, 0.25),
#             weighted_quantile(vals, w, 0.75)
#         )
#     # Sort by rel_time
#     rel_times = sorted(stats.keys())
#     means = [stats[t][0] for t in rel_times]
#     q25s = [stats[t][1] for t in rel_times]
#     q75s = [stats[t][2] for t in rel_times]
#     return rel_times, means, q25s, q75s

def calculate_weighted_stats(df, varname):
    """Returns dict: rel_time -> (mean, q25, q75) for the given variable, unweighted."""
    stats = {}
    for t, group in df.groupby('rel_time'):
        vals = group[varname]
        stats[t] = (
            vals.mean(),
            vals.quantile(0.25),
            vals.quantile(0.75)
        )
    # Sort by rel_time
    rel_times = sorted(stats.keys())
    means = [stats[t][0] for t in rel_times]
    q25s = [stats[t][1] for t in rel_times]
    q75s = [stats[t][2] for t in rel_times]
    return rel_times, means, q25s, q75s

def plot_weighted_stats(stats_dict, th_limits, region, airstream, varname):
    # plt.figure(figsize=(8, 6)) #old: (13, 8)
    
    # test
    #####
    pt = 1./72.27 # Hundreds of years of history... 72.27 points to an inch.

    textwidth = 418.25368


    figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
                # CQG is only one column
                # Add more journals below. Can add more properties to each journal

    my_width = figsizes["10"]
    # Our figure's aspect ratio
    golden = (1 + 5 ** 0.5) / 2
    fourto3 = 4/3

    plt.figure(figsize = (my_width, my_width/fourto3))
    #####

    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        plt.plot(rel_times, means, label=f'mean {varname} ({timeframe})', color=colors[timeframe])
        plt.fill_between(rel_times, q25s, q75s, alpha=0.2, color=colors[timeframe], label=f'IQR ({timeframe})')
    plt.xlabel('hours relative to CAO [h]')
    if varname == 'Q':
        plt.ylabel(f'q [g kg${-1}$]')
    elif varname == 'TH':
        plt.ylabel(f'$\Theta$ [K]')
    else:
        plt.ylabel(f'{varname} [{units_dict[varname]}]')
    plt.xlim(-168, 0)
    plt.xticks(range(-168, 1, 12))
    #### new
    plt.gca().set_xticks(np.arange(-168, 1, 24))
    plt.gca().set_xticklabels([str(int(x//24)) for x in np.arange(-168, 1, 24)])
    plt.gca().set_xticks(np.arange(-168, 1, 6), minor=True)
    ####
    plt.xlabel('days relative to CAO')
    plt.yticks()
    plt.ylim(th_limits[region])
    if varname == 'TH':
        plt.ylim(250, 300)
    plt.grid(True, alpha=0.3)
    plt.legend()
    # plt.title(f'{airstream} - {region}: Mean and IQR of {varname} vs. Relative Time')
    plt.tight_layout()
    plt.savefig(f"{figpath}evo_{varname}_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
    plt.show()

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
    
# %%
def plot_combined_evolution_rates(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    """
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"{config.LAGRANTO_DIR}/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
        dfs[timeframe] = pd.read_csv(path)
    
    # Calculate evolution stats
    stats_dict = {}
    for timeframe in timeframes:
        df = dfs[timeframe]
        df_air = df[df['airstream'] == airstream]
        if df_air.empty:
            continue
        stats_dict[timeframe] = calculate_weighted_stats(df_air, varname)
    
    # Prepare rate data
    df_air_present = df_present[df_present['airstream'] == airstream]
    df_air_future = df_future[df_future['airstream'] == airstream]
    
    if df_air_present.empty or df_air_future.empty:
        print(f"No data for {airstream} in {region}")
        return
    
    # Set up the figure with custom width ratios
    pt = 1./72.27
    textwidth = 418.25368
    fig_width = textwidth * pt  # Use textwidth for overall figure width
    fourto3 = 4/3
    figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
    my_width = figsizes["10"]
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    # fig, axes = plt.subplots(1, 3, figsize=((my_width*3, my_width/fourto3)), 
    #                         gridspec_kw={'width_ratios': [3, 1, 1]})
    fig, axes = plt.subplots(1, 5, figsize=((my_width/3*7, my_width/fourto3)), 
                            gridspec_kw={'width_ratios': [3, 1, 1, 1, 1]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First subplot: Evolution plot
    ax1 = axes[0]
    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})', color=colors[timeframe])
        ax1.fill_between(rel_times, q25s, q75s, alpha=0.2, color=colors[timeframe], label=f'IQR ({timeframe})')
        # Always use TH to find the minimum time, regardless of varname
        if varname != 'TH':
            # Need to calculate TH stats to find minimum TH time
            th_stats = calculate_weighted_stats(dfs[timeframe][dfs[timeframe]['airstream'] == airstream], 'TH')
            th_rel_times, th_means, _, _ = th_stats
            min_th_idx = np.argmin(th_means)
            min_th_time = th_rel_times[min_th_idx]
        else:
            # For TH, use the current stats
            min_th_idx = np.argmin(means)
            min_th_time = rel_times[min_th_idx]
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--', alpha=0.7)
    
    ax1.set_xlabel('days relative to CAO')
    if varname == 'Q':
        ax1.set_ylabel(f'q [g kg$^{{-1}}$]')
        ax1.set_ylim(0, 5)
    elif varname == 'TH':
        ax1.set_ylabel(f'$\\Theta$ [K]')
        ax1.set_ylim(250, 300)
    else:
        ax1.set_ylabel(f'{varname} [{units_dict[varname]}]')
    
    ax1.set_xlim(-168, 0)
    ax1.set_xticks(np.arange(-168, 1, 24))
    ax1.set_xticklabels([str(int(x//24)) for x in np.arange(-168, 1, 24)])
    ax1.set_xticks(np.arange(-168, 1, 6), minor=True)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Second subplot: Cooling rates boxplot
    ax2 = axes[1]
    if varname == 'TH':
        cooling_rate_present = df_air_present['cooling_rate_6h'].dropna()
        cooling_rate_future = df_air_future['cooling_rate_6h'].dropna()
        rate_label = '$\\Theta$ rate [K 6h$^{-1}$]'
    elif varname == 'Q':
        cooling_rate_present = df_air_present['q_rate_6h_cooling'].dropna()
        cooling_rate_future = df_air_future['q_rate_6h_cooling'].dropna()
        rate_label = 'q rate [g kg$^{-1}$ 6h$^{-1}$]'

    data_cooling = [cooling_rate_present, cooling_rate_future]

    # Vertical boxplots
    box1 = ax2.boxplot([data_cooling[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=1.5),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))

    box2 = ax2.boxplot([data_cooling[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=1.5),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels(['present', 'future'])
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)
    ax2.set_title('cooling phase')

    # Third subplot: WBC Warming rates boxplot  
    ax3 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_wbc'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_wbc'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming_wbc'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming_wbc'].dropna()

    data_warming = [warming_rate_present, warming_rate_future]

    box3 = ax3.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=1.5),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))

    box4 = ax3.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=1.5),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))

    ax3.set_xticks([1, 2])
    ax3.set_xticklabels(['present', 'future'])
    # Don't set ylabel for ax3 - it will share the same label as ax2
    ax3.grid(True, axis='y', alpha=0.3)
    ax3.set_title('warming phase - WBC')

    # Fourth Subplot: non-WBC Warming rates boxplot
    ax4 = axes[3]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_non_wbc'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_non_wbc'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming_non_wbc'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming_non_wbc'].dropna()
    
    data_warming = [warming_rate_present, warming_rate_future]

    box5 = ax4.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=1.5),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))
    
    box6 = ax4.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=1.5),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))
    
    ax4.set_xticks([1, 2])
    ax4.set_xticklabels(['present', 'future'])
    # Don't set ylabel for ax4 - it will share the same label as ax2
    ax4.grid(True, axis='y', alpha=0.3)
    ax4.set_title('warming phase - non-WBC')

    # Fifth subplot: Warming rates boxplot
    ax5 = axes[4]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=1.5),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=1.5),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=1.5))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels(['present', 'future'])
    ax5.grid(True, axis='y', alpha=0.3)
    ax5.set_title('warming phase - all trajectories')

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
    plt.show()

def plot_th_min_map(region, airstream, timeframe):
    """
    Plot the geographic locations where TH minimum is reached on a map
    """
    # Load the summary data
    path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
    df = pd.read_csv(path)
    
    # Filter for specific airstream
    df_air = df[df['airstream'] == airstream].copy()
    
    if df_air.empty:
        print(f"No data for {airstream} in {region} ({timeframe})")
        return
    
    # Remove NaN values
    df_clean = df_air.dropna(subset=['TH_min_lon', 'TH_min_lat'])
    
    print(f"Plotting {len(df_clean)} TH minimum locations for {airstream} in {region} ({timeframe})")
    
    # Set up the map projection and extent based on region
    if region == 'NA':
        # North Atlantic
        projection = ccrs.PlateCarree()
        extent = [-80, 20, 35, 75]  # [lon_min, lon_max, lat_min, lat_max]
    else:  # JS
        # Japan Sea
        projection = ccrs.PlateCarree()
        extent = [120, 150, 35, 55]
    
    # Create the map
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    #ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Add map features
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, color='lightgray')


    ax.add_feature(cfeature.OCEAN, color='lightblue')
    
    # Add gridlines
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    
    # Scatter plot of TH minimum locations
    scatter = ax.scatter(df_clean['TH_min_lon'], df_clean['TH_min_lat'], 
                        c=df_clean['TH_min'], cmap='coolwarm', 
                        s=30, alpha=0.7, edgecolors='black', linewidth=0.3,
                        transform=ccrs.PlateCarree())
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('TH minimum [K]', fontsize=12)
    
    # Add title
    plt.title(f'TH Minimum Locations: {airstream} - {region} ({timeframe})\n'
             f'N = {len(df_clean)} trajectories', fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(f"{figpath}th_min_map_{airstream}_{region}_{timeframe}.png", 
                dpi=300, bbox_inches='tight')
    plt.show()

# %% NEW
# For each traj_id, find TH at rel_time=0, TH at rel_time=-168, and minimum TH
regions = ['JS', 'NA']
timeframes = ['present', 'future']
airstreams = ['LL', 'UL']

regions = ['NA']  # For testing, can be expanded later
timeframes = ['present']
airstreams = ['LL']

for region in regions:
    for timeframe in timeframes:
        path = f"{config.LAGRANTO_DIR}/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
        df = pd.read_csv(path)
        summary_list = []
        for airstream in airstreams:
            df_air = df[df['airstream'] == airstream]
            if df_air.empty:
                continue
            grouped = df_air.groupby('traj_id')
            for traj_id, group in grouped:
                # Sort by rel_time to ensure proper ordering
                group = group.sort_values('rel_time')
                
                th_0 = group.loc[np.isclose(group['rel_time'], 0.0), 'TH']
                th_168 = group.loc[np.isclose(group['rel_time'], -168.0), 'TH']
                th_min = group['TH'].min()
                TH_min_time = group.loc[group['TH'] == th_min, 'rel_time'].values[0] if not group[group['TH'] == th_min].empty else np.nan
                TH_min_lon = group.loc[group['TH'] == th_min, 'lon'].values[0] if not group[group['TH'] == th_min].empty else np.nan
                TH_min_lat = group.loc[group['TH'] == th_min, 'lat'].values[0] if not group[group['TH'] == th_min].empty else np.nan
                
                # Calculate 6-hourly change rates for each time step
                group = group.copy()  # Avoid SettingWithCopyWarning
                group['TH_6h_change_rate'] = np.nan
                group['Q_6h_change_rate'] = np.nan
                for i in range(1, len(group)):
                    dt = group.iloc[i]['rel_time'] - group.iloc[i-1]['rel_time']
                    if dt == 6:  # Only calculate for 6-hour intervals
                        dTH = group.iloc[i]['TH'] - group.iloc[i-1]['TH']
                        group.iloc[i, group.columns.get_loc('TH_6h_change_rate')] = dTH  # Change per 6 hours
                        dQ = group.iloc[i]['Q'] - group.iloc[i-1]['Q']
                        group.iloc[i, group.columns.get_loc('Q_6h_change_rate')] = dQ  # Change per 6 hours

                # Find cooling and warming periods based on TH_min_time
                # Cooling: from -168h to TH_min_time
                cooling_mask = (group['rel_time'] >= -168) & (group['rel_time'] <= TH_min_time) & ~group['TH_6h_change_rate'].isna()
                cooling_rates = group.loc[cooling_mask, 'TH_6h_change_rate']
                q_rates_cooling = group.loc[cooling_mask, 'Q_6h_change_rate']
                cooling_rate_6h = cooling_rates.mean() if len(cooling_rates) > 0 else np.nan
                q_rate_6h_cooling = q_rates_cooling.mean() if len(q_rates_cooling) > 0 else np.nan
                
                # Warming: from TH_min_time to 0h, only consider SST > 0 (over open water)
                warming_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna()
                sst_positive_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna() & (group['SST'] > 0)
                warming_rates = group.loc[warming_mask, 'TH_6h_change_rate']
                q_rates_warming = group.loc[warming_mask, 'Q_6h_change_rate']
                warming_rate_6h = warming_rates.mean() if len(warming_rates) > 0 else np.nan
                q_rate_6h_warming = q_rates_warming.mean() if len(q_rates_warming) > 0 else np.nan

                #Warming: WBC vs non-WBC
                if region == 'NA':
                    non_WBC_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna() & (group['SST'] > 0) & (group['lat'] > 50)
                    WBC_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna() & (group['SST'] > 0) & (group['lat'] <= 50)
                if region == 'JS':
                    non_WBC_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna() & (group['SST'] > 0) & (group['lat'] > 65)
                    WBC_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna() & (group['SST'] > 0) & (group['lat'] <= 65)
                warming_rates_non_wbc = group.loc[non_WBC_mask, 'TH_6h_change_rate']
                warming_rate_6h_non_wbc = warming_rates_non_wbc.mean() if not warming_rates_non_wbc.empty else np.nan
                q_rates_non_wbc = group.loc[non_WBC_mask, 'Q_6h_change_rate']
                q_rate_6h_non_wbc = q_rates_non_wbc.mean() if not q_rates_non_wbc.empty else np.nan
                warming_rates_wbc = group.loc[WBC_mask, 'TH_6h_change_rate']
                warming_rate_6h_wbc = warming_rates_wbc.mean() if not warming_rates_wbc.empty else np.nan
                q_rates_wbc = group.loc[WBC_mask, 'Q_6h_change_rate']
                q_rate_6h_wbc = q_rates_wbc.mean() if not q_rates_wbc.empty else np.nan

                summary_list.append({
                    'region': region,
                    'timeframe': timeframe,
                    'airstream': airstream,
                    'traj_id': traj_id,
                    'TH_at_0': th_0.iloc[0] if not th_0.empty else np.nan,
                    'TH_at_-168': th_168.iloc[0] if not th_168.empty else np.nan,
                    'TH_min': th_min,
                    'TH_min_time': TH_min_time,
                    'TH_min_lon': TH_min_lon,
                    'TH_min_lat': TH_min_lat,
                    'cooling_rate_6h': cooling_rate_6h,  # NEW: Average K per 6-hour interval during cooling
                    'warming_rate_6h': warming_rate_6h,   # NEW: Average K per 6-hour interval during warming
                    'warming_rate_6h_non_wbc': warming_rate_6h_non_wbc,  # NEW: Average K per 6-hour interval during warming (non-WBC)
                    'warming_rate_6h_wbc': warming_rate_6h_wbc,   # NEW: Average K per 6-hour interval during warming (WBC)
                    'q_rate_6h_cooling': q_rate_6h_cooling,  # NEW: Average Q per 6-hour interval during cooling
                    'q_rate_6h_warming': q_rate_6h_warming,   # NEW: Average Q per 6-hour interval during warming
                    'q_rate_6h_non_wbc': q_rate_6h_non_wbc,  # NEW: Average Q per 6-hour interval during warming (non-WBC)
                    'q_rate_6h_wbc': q_rate_6h_wbc    # NEW: Average Q per 6-hour interval during warming (WBC)
                })

        summary_df = pd.DataFrame(summary_list)

        # for each traj_id, get the slope of cooling and warming
        summary_df['cooling'] = summary_df['TH_min'] - summary_df['TH_at_-168']
        summary_df['warming'] = summary_df['TH_at_0'] - summary_df['TH_min']
        summary_df['cooling_time'] = summary_df['TH_min_time'] + 168
        summary_df['warming_time'] = 0-summary_df['TH_min_time']

        # Example: save to CSV
        summary_df.to_csv(f"{datapath}TH_summary_{region}_{timeframe}.csv", index=False, header=True)
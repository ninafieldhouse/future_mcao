# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm


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

# plt.rcParams["font.family"] = "sans-serif"
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
airstreams = ['BCAR', 'DI', 'nonBCAR'
              , 'LL'
              ]
airstreams = ['LL', 'UL']
timeframes = ['present', 'future']
colors = {'present': 'blue', 'future': 'orange'}
figpath = '/home/nfieldhouse/ma/plots/traj/'
datapath = '/net/litho/atmosdyn2/nfieldhouse/data/lagranto/'

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
    t_stat : float
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
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    ax2.set_xticklabels([f'present\n{data_cooling[0].mean():.2f}', f'future\n{data_cooling[1].mean():.2f}'])
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)
    ax2.set_title('cooling phase')

    # Third subplot: WBC Warming rates boxplot  
    ax3 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_wbc'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_wbc'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_wbc'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_wbc'].dropna()

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
    ax3.set_xticklabels([f'present\n{data_warming[0].mean():.2f}', f'future\n{data_warming[1].mean():.2f}'])
    # Don't set ylabel for ax3 - it will share the same label as ax2
    ax3.grid(True, axis='y', alpha=0.3)
    ax3.set_title('warming phase - WBC')

    # Fourth Subplot: non-WBC Warming rates boxplot
    ax4 = axes[3]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_non_wbc'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_non_wbc'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_non_wbc'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_non_wbc'].dropna()
    
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
    ax4.set_xticklabels([f'present\n{data_warming[0].mean():.2f}', f'future\n{data_warming[1].mean():.2f}'])
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
    ax5.set_xticklabels([f'present\n{data_warming[0].mean():.2f}', f'future\n{data_warming[1].mean():.2f}'])
    ax5.grid(True, axis='y', alpha=0.3)
    ax5.set_title('warming phase - all trajectories')

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}.pdf", bbox_inches='tight')
    plt.show()
def plot_combined_evolution_rates_cooling_warming(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    """

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4      # Y-axis tick font size
})
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    my_width = 6.69291
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    # fig, axes = plt.subplots(1, 3, figsize=((my_width*3, my_width/fourto3)), 
    #                         gridspec_kw={'width_ratios': [3, 1, 1]})
    fig, axes = plt.subplots(1, 3, figsize=(my_width, 2.2), 
                            gridspec_kw={'width_ratios': [3, 1, 1,]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First subplot: Evolution plot
    ax1 = axes[0]
    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})' if varname != 'TH' else f'mean $\\Theta$ ({timeframe})', color=colors[timeframe], linewidth=0.8)
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
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--', alpha=0.7, linewidth=0.8)
    
    ax1.set_xlabel('days relative to t$_{MCAO}$')
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
    # ax1.set_title(f'a) {region} - {airstream}', loc='left')
    
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

    is_sig_cooling, p_value_cooling, t_stat_cooling = compute_significance_ttest(cooling_rate_present.values, cooling_rate_future.values)
    data_cooling = [cooling_rate_present, cooling_rate_future]
    print(f"Cooling rates - {region} - {airstream} - {varname}: p-value = {p_value_cooling}, t-stat = {t_stat_cooling}")

    asterisk_cooling = '*' if is_sig_cooling else ''

    # Vertical boxplots
    box1 = ax2.boxplot([data_cooling[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_cooling[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels([f'present\n {data_cooling[0].mean():.2f}', f'future{asterisk_cooling}\n {data_cooling[1].mean():.2f}'], fontsize=6)
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)
    # ax2.set_title('b) cooling phase', loc='left')

    # Fifth subplot: Warming rates boxplot
    ax5 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    print(f"Warming rates - {region} - {airstream} - {varname}: p-value = {p_value_warming}, t-stat = {t_stat_warming}")
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels([f'present\n {data_warming[0].mean():.2f}', f'future{asterisk_warming}\n {data_warming[1].mean():.2f}'], fontsize=6)
    ax5.grid(True, axis='y', alpha=0.3)
    # ax5.set_title('c) warming phase', loc='left')

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming.pdf", bbox_inches='tight')
    plt.show()
def plot_combined_evolution_rates_wbc_nonwbc_land(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Vertical boxplot for wbc warming
    2. Vertical boxplot for non-wbc warming
    3. Vertical boxplot for land warming
    """

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4      # Y-axis tick font size
})
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    my_width = 6.69291

    if region == 'JS':
        letters = ['a', 'b', 'c']
    elif region == 'NA':
        letters = ['d', 'e', 'f']
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    # fig, axes = plt.subplots(1, 3, figsize=((my_width*3, my_width/fourto3)), 
    #                         gridspec_kw={'width_ratios': [3, 1, 1]})
    fig, axes = plt.subplots(1, 3, figsize=(my_width/5*3, 2.2), 
                            gridspec_kw={'width_ratios': [1, 1, 1,]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First Subplot: WBC warming rates boxplot
    ax1 = axes[0]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_wbc'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_wbc'].dropna()
        rate_label = '$\\Theta$ rate [K 6h$^{-1}$]'
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_wbc'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_wbc'].dropna()
        rate_label = 'q rate [g kg$^{-1}$ 6h$^{-1}$]'
   
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    print(f"WBC Warming rates - {region} - {airstream} - {varname}: p-value = {p_value_warming}, t-stat = {t_stat_warming}")
    print(f'WBC: number of present trajectories: {len(warming_rate_present)}, number of future trajectories: {len(warming_rate_future)}')
    box1 = ax1.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box2 = ax1.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax1.set_xticks([1, 2])
    ax1.set_xticklabels([f'present\n{data_warming[0].mean():.2f}\nN = {len(warming_rate_present)}', f'future{asterisk_warming}\n{data_warming[1].mean():.2f}\nN = {len(warming_rate_future)}'], fontsize=6)
    ax1.set_ylabel(rate_label)
    ax1.grid(True, axis='y', alpha=0.3)
    ax1.set_title(f'{letters[0]}) {region} - {airstream} - WBC', loc='left')


    
    # Second subplot: Non-WBC warming rates boxplot
    ax2 = axes[1]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_non_wbc'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_non_wbc'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_non_wbc'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_non_wbc'].dropna()

    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    print(f"Non-WBC Warming rates - {region} - {airstream} - {varname}: p-value = {p_value_warming}, t-stat = {t_stat_warming}")
    print(f'Non-WBC: number of present trajectories: {len(warming_rate_present)}, number of future trajectories: {len(warming_rate_future)}')
    # Vertical boxplots
    box1 = ax2.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels([f'present\n{data_warming[0].mean():.2f}\nN = {len(warming_rate_present)}', f'future{asterisk_warming}\n{data_warming[1].mean():.2f}\nN = {len(warming_rate_future)}'], fontsize=6)
    ax2.grid(True, axis='y', alpha=0.3)
    ax2.set_title(f'{letters[1]}) {region} - {airstream} - non-WBC', loc='left')

    # Third subplot: Warming rates land boxplot
    ax5 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h_land'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h_land'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_land'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_land'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    print(f"Land Warming rates - {region} - {airstream} - {varname}: p-value = {p_value_warming}, t-stat = {t_stat_warming}")
    print(f'Land Warming rates - {region} - {airstream} - {varname}: number of present trajectories: {len(warming_rate_present)}, number of future trajectories: {len(warming_rate_future)}')
    
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels([f'present\n{data_warming[0].mean():.2f}\nN = {len(warming_rate_present)}', f'future{asterisk_warming}\n{data_warming[1].mean():.2f}\nN = {len(warming_rate_future)}'], fontsize=6)
    ax5.grid(True, axis='y', alpha=0.3)
    ax5.set_title(f'{letters[2]}) {region} - {airstream} - land', loc='left')

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_warmingphases.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_warmingphases.pdf", bbox_inches='tight')
    plt.show()

def plot_combined_evolution_rates_cooling_warming_min24h(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    
    Filters out cooling times shorter than 24 hours.
    """

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4      # Y-axis tick font size
})
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    
    # Filter out cooling times shorter than 24 hours
    df_air_present = df_air_present[df_air_present['cooling_time'] >= 24]
    df_air_future = df_air_future[df_air_future['cooling_time'] >= 24]
    
    # Set up the figure with custom width ratios
    pt = 1./72.27
    textwidth = 418.25368
    fig_width = textwidth * pt  # Use textwidth for overall figure width
    fourto3 = 4/3
    figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
    my_width = 6.69291
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    fig, axes = plt.subplots(1, 3, figsize=(my_width, 2.2), 
                            gridspec_kw={'width_ratios': [3, 1, 1,]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First subplot: Evolution plot
    ax1 = axes[0]
    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})', color=colors[timeframe], linewidth=0.8)
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
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--', alpha=0.7, linewidth=0.8)
    
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

    is_sig_cooling, p_value_cooling, t_stat_cooling = compute_significance_ttest(cooling_rate_present.values, cooling_rate_future.values)
    data_cooling = [cooling_rate_present, cooling_rate_future]

    asterisk_cooling = '*' if is_sig_cooling else ''

    # Vertical boxplots
    box1 = ax2.boxplot([data_cooling[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_cooling[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels([f'present\n mean = {data_cooling[0].mean():.2f}\n median = {data_cooling[0].median():.2f}', f'future\n mean{asterisk_cooling} = {data_cooling[1].mean():.2f}\n median = {data_cooling[1].median():.2f}'], fontsize=6)
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)

    # Fifth subplot: Warming rates boxplot
    ax5 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels([f'present\n mean = {data_warming[0].mean():.2f}\n median = {data_warming[0].median():.2f}', f'future\n mean{asterisk_warming} = {data_warming[1].mean():.2f}\n median = {data_warming[1].median():.2f}'], fontsize=6)
    ax5.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_min24h.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_min24h.pdf", bbox_inches='tight')
    plt.show()

def plot_combined_evolution_rates_cooling_warming_no_outliers(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    
    Filters out the top and bottom 1% of cooling rates.
    """

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4      # Y-axis tick font size
})
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    
    # Filter out top and bottom 1% of cooling rates
    if varname == 'TH':
        cooling_col = 'cooling_rate_6h'
    elif varname == 'Q':
        cooling_col = 'q_rate_6h_cooling'
    
    p1_present = df_air_present[cooling_col].quantile(0.01)
    p99_present = df_air_present[cooling_col].quantile(0.99)
    p1_future = df_air_future[cooling_col].quantile(0.01)
    p99_future = df_air_future[cooling_col].quantile(0.99)
    
    df_air_present = df_air_present[(df_air_present[cooling_col] >= p1_present) & (df_air_present[cooling_col] <= p99_present)]
    df_air_future = df_air_future[(df_air_future[cooling_col] >= p1_future) & (df_air_future[cooling_col] <= p99_future)]
    
    # Set up the figure with custom width ratios
    pt = 1./72.27
    textwidth = 418.25368
    fig_width = textwidth * pt  # Use textwidth for overall figure width
    fourto3 = 4/3
    figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
    my_width = 6.69291
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    fig, axes = plt.subplots(1, 3, figsize=(my_width, 2.2), 
                            gridspec_kw={'width_ratios': [3, 1, 1,]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First subplot: Evolution plot
    ax1 = axes[0]
    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})', color=colors[timeframe], linewidth=0.8)
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
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--', alpha=0.7, linewidth=0.8)
    
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

    is_sig_cooling, p_value_cooling, t_stat_cooling = compute_significance_ttest(cooling_rate_present.values, cooling_rate_future.values)
    data_cooling = [cooling_rate_present, cooling_rate_future]

    asterisk_cooling = '*' if is_sig_cooling else ''

    # Vertical boxplots
    box1 = ax2.boxplot([data_cooling[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_cooling[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels([f'present\n mean = {data_cooling[0].mean():.2f}\n median = {data_cooling[0].median():.2f}', f'future\n mean{asterisk_cooling} = {data_cooling[1].mean():.2f}\n median = {data_cooling[1].median():.2f}'], fontsize=6)
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)

    # Fifth subplot: Warming rates boxplot
    ax5 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels([f'present\n mean = {data_warming[0].mean():.2f}\n median = {data_warming[0].median():.2f}', f'future\n mean{asterisk_warming} = {data_warming[1].mean():.2f}\n median = {data_warming[1].median():.2f}'], fontsize=6)
    ax5.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_no_outliers.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_no_outliers.pdf", bbox_inches='tight')
    plt.show()

def plot_th_min_map(region, airstream, timeframe, only_non_WBC=False, colorby='TH_min_time'):
    """
    Plot the geographic locations where TH minimum is reached on a map
    """
    # Load the summary data
    path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
    df = pd.read_csv(path)
    
    # Filter for specific airstream
    df_air = df[df['airstream'] == airstream].copy()
    if colorby == 'TH_min':
        if region == 'NA':
            cbar_min = 250
            cbar_max = 270
        elif region == 'JS':
            cbar_min = 250  
            cbar_max = 260
        # cbar_min = df_air['TH_min'].quantile(0.25)
        # cbar_max = df_air['TH_min'].quantile(0.75)
    else:
        cbar_min = -168
        cbar_max = 0
    
    if df_air.empty:
        print(f"No data for {airstream} in {region} ({timeframe})")
        return

    if only_non_WBC:
        # select rows where 'warming_rate_6h_non_WBC' is not NaN
        df_air = df_air[~df_air['warming_rate_6h_non_wbc'].isna()]
    
    # Filter for IQR of TH minimum
    th_min_q25 = df_air['TH_min'].quantile(0.25)
    th_min_q75 = df_air['TH_min'].quantile(0.75)
    df_air = df_air[(df_air['TH_min'] >= th_min_q25) & (df_air['TH_min'] <= th_min_q75)]
    
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

    projection = ccrs.Orthographic(central_longitude=0, central_latitude=90)
    extent = [-180, 180, 20, 90]

    
    # Create the map
    fig = plt.figure(figsize=(10, 10))
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

    if colorby == 'TH_min_time':
        cmap = 'viridis'
        levels = np.linspace(cbar_min, cbar_max, 15)
        norm = BoundaryNorm(boundaries=levels, ncolors=256)
        vmin = None
        vmax = None
    elif colorby == 'TH_min':
        cmap = 'RdBu_r'
        norm = None
        vmin = cbar_min
        vmax = cbar_max
    
    
    # Scatter plot of TH minimum locations
    scatter = ax.scatter(df_clean['TH_min_lon'], df_clean['TH_min_lat'], 
                        c=df_clean[colorby], cmap=cmap, vmin=vmin, vmax=vmax, norm=norm,
                        s=30, alpha=0.7, edgecolors='black', linewidth=0.3,
                        transform=ccrs.PlateCarree())
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('TH minimum [K]', fontsize=12)
    
    # Add title
    plt.title(f'TH Minimum Locations: {airstream} - {region} ({timeframe})\n'
             f'N = {len(df_clean)} trajectories (IQR), only non-WBC = {only_non_WBC}', fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(f"{figpath}th_min_map_{airstream}_{region}_{timeframe}_onlyNonWBC_{only_non_WBC}.png", 
                dpi=300, bbox_inches='tight')
    plt.show()

# %%
def plot_th_min_map_grid(region, colorby='TH_min_time'):
    """
    Create a 2x2 grid of TH minimum location maps
    Rows: present (top), future (bottom)  
    Columns: all trajectories (left), only non-WBC (right)
    """
    # Set up colorbar limits
    if colorby == 'TH_min':
        if region == 'NA':
            cbar_min = 250
            cbar_max = 270
        elif region == 'JS':
            cbar_min = 250  
            cbar_max = 260
        cmap = 'RdBu_r'
        cbar_label = 'TH minimum [K]'
    else:  # TH_min_time
        cbar_min = -168
        cbar_max = 0
        cmap = 'viridis'
        cbar_label = 'TH minimum time [h]'
    
    
    projection = ccrs.Orthographic(central_longitude=0, central_latitude=90)
    extent = [-180, 180, 15, 90]
    
    # Create figure and subplots
    fig = plt.figure(figsize=(16, 16))
    
    # Define subplot parameters
    timeframes = ['present', 'future']
    only_non_WBC_vals = [False, True]
    airstream = 'LL'  # You can make this a parameter if needed
    
    axes = []
    for i, timeframe in enumerate(timeframes):
        for j, only_non_WBC in enumerate(only_non_WBC_vals):
            ax = fig.add_subplot(2, 2, i*2 + j + 1, projection=projection)
            axes.append(ax)
            
            # Load and filter data
            path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
            df = pd.read_csv(path)
            df_air = df[df['airstream'] == airstream].copy()
            
            if df_air.empty:
                print(f"No data for {airstream} in {region} ({timeframe})")
                continue
            
            if only_non_WBC:
                # Select rows where 'warming_rate_6h_non_wbc' is not NaN
                df_air = df_air[~df_air['warming_rate_6h_non_wbc'].isna()]
            
            # Filter for IQR of TH minimum
            th_min_q25 = df_air['TH_min'].quantile(0.25)
            th_min_q75 = df_air['TH_min'].quantile(0.75)
            df_air = df_air[(df_air['TH_min'] >= th_min_q25) & (df_air['TH_min'] <= th_min_q75)]
            
            # Remove NaN values
            df_clean = df_air.dropna(subset=['TH_min_lon', 'TH_min_lat'])
                        
            # Add map features
            ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, color='lightgray')
            ax.add_feature(cfeature.OCEAN, color='lightblue')
            
            # Add gridlines
            gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
            if region == 'JS':
                gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                                 xlocs=range(-180, 181, 20), ylocs=range(30, 81, 10))
            gl.top_labels = False
            gl.right_labels = False
            # Only show labels on left and bottom edges
            if j > 0:  # Right column
                gl.left_labels = False
            if i == 0:  # Top row
                gl.bottom_labels = False
            
            # Set up colormap and normalization
            if colorby == 'TH_min_time':
                levels = np.linspace(cbar_min, cbar_max, 15)
                norm = BoundaryNorm(boundaries=levels, ncolors=256)
                vmin = None
                vmax = None
            else:  # TH_min
                norm = None
                vmin = cbar_min
                vmax = cbar_max
            
            # Scatter plot of TH minimum locations
            scatter = ax.scatter(df_clean['TH_min_lon'], df_clean['TH_min_lat'], 
                               c=df_clean[colorby], cmap=cmap, vmin=vmin, vmax=vmax, norm=norm,
                               s=30, alpha=0.7, edgecolors='black', linewidth=0.3,
                               transform=ccrs.PlateCarree())

            ax.set_extent(extent, crs=ccrs.PlateCarree())
            
            # Add subplot titles (removed fontsize)
            wbc_label = 'only non-WBC' if only_non_WBC else 'all trajectories'
            ax.set_title(f'{timeframe} - {wbc_label}\nN = {len(df_clean)}', pad=10)
    
    # Create shared colorbar at bottom
    # Adjust subplot spacing to make room for colorbar
    plt.subplots_adjust(bottom=0.15, hspace=0.3, wspace=0.1)
    
    # Add colorbar (removed fontsize)
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.03])  # [left, bottom, width, height]
    cbar = fig.colorbar(scatter, cax=cbar_ax, orientation='horizontal')
    cbar.set_label(cbar_label)
    
    # Add main title (removed fontsize and y parameter)
    fig.suptitle(f'TH Minimum Locations: {airstream} - {region}')
    
    # Save and show
    plt.savefig(f"{figpath}th_min_map_grid_{airstream}_{region}_{colorby}.png", 
                dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}th_min_map_grid_{airstream}_{region}_{colorby}.pdf", 
                bbox_inches='tight')  # PDF format for vector graphics
    plt.show()
    
    return fig, axes

# %%
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
def plot_th_min_map_dbscan(region, airstream, timeframe, only_non_WBC = False,eps=2.0, min_samples=5):
    """
    Plot the geographic locations where TH minimum is reached on a map with DBSCAN clustering
    """
    from sklearn.cluster import DBSCAN
    
    # Load the summary data
    path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
    df = pd.read_csv(path)
    
    # Filter for specific airstream
    df_air = df[df['airstream'] == airstream].copy()
    
    if df_air.empty:
        print(f"No data for {airstream} in {region} ({timeframe})")
        return
    
    if only_non_WBC:
        # select rows where 'warming_rate_6h_non_WBC' is not NaN
        df_air = df_air[~df_air['warming_rate_6h_non_wbc'].isna()]
    
    # Filter for IQR of TH minimum
    th_min_q25 = df_air['TH_min'].quantile(0.25)
    th_min_q75 = df_air['TH_min'].quantile(0.75)
    df_air = df_air[(df_air['TH_min'] >= th_min_q25) & (df_air['TH_min'] <= th_min_q75)]
    
    # Remove NaN values
    df_clean = df_air.dropna(subset=['TH_min_lon', 'TH_min_lat'])
    
    print(f"Plotting {len(df_clean)} TH minimum locations for {airstream} in {region} ({timeframe})")
    
    # Perform DBSCAN clustering
    coordinates = df_clean[['TH_min_lon', 'TH_min_lat']].values
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    cluster_labels = dbscan.fit_predict(coordinates)
    df_clean['cluster'] = cluster_labels
    
    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = list(cluster_labels).count(-1)
    
    # Calculate extent based on cluster points only (excluding noise)
    cluster_points = df_clean[df_clean['cluster'] != -1]
    
    # Set up the map projection and extent based on region
    if region == 'NA':
        extent = [-180, 0, 30, 90]
        projection = ccrs.PlateCarree()
    else:  # JS - spanning across dateline
        # For JS: from 0° to 140°W (which is -140°)
        # This spans across the dateline: 0° to 180°, then -180° to -140°
        projection = ccrs.PlateCarree(central_longitude = 110)
        extent = [0, 220, 30, 90]  # [lon_min, lon_max, lat_min, lat_max]

    projection = ccrs.Orthographic(central_longitude=0, central_latitude=90)
    
    # Create the map with larger figure size
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    
    # Set the extent to focus on the region
    extent = [-180, 180, 45, 90]
   # ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    # Add map features
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, color='lightgray')
    ax.add_feature(cfeature.OCEAN, color='lightblue')
    
    # Add gridlines with special handling for dateline-spanning maps
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                         xlocs=range(-180, 181, 20), ylocs=range(30, 81, 10))
    
    gl.top_labels = False
    gl.right_labels = False
    
    # Create colormap for clusters
    import matplotlib.cm as cm
    if n_clusters > 0:
        cluster_colors = cm.viridis(np.linspace(0, 1, n_clusters))
    
    # Plot clusters
    unique_labels = set(cluster_labels)
    for k in unique_labels:
        if k == -1:
            # Noise points in black
            cluster_data = df_clean[df_clean['cluster'] == k]
            if len(cluster_data) > 0:
                ax.scatter(cluster_data['TH_min_lon'], cluster_data['TH_min_lat'], 
                          c='black', s=20, alpha=0.5, 
                          transform=ccrs.PlateCarree(), label='Noise')
        else:
            # Cluster points
            cluster_data = df_clean[df_clean['cluster'] == k]
            if len(cluster_data) > 0:
                ax.scatter(cluster_data['TH_min_lon'], cluster_data['TH_min_lat'], 
                          c=[cluster_colors[k]], s=30, alpha=0.7, 
                          edgecolors='black', linewidth=0.3,
                          transform=ccrs.PlateCarree(), label=f'Cluster {k+1}')
    
    # Add legend with better positioning
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    
    # Add title with adjusted padding
    plt.title(f'TH Minimum Locations (DBSCAN): {airstream} - {region} ({timeframe})\n'
             f'N = {len(df_clean)} trajectories (IQR), {n_clusters} clusters, {n_noise} noise points',pad=15)
    
    # Use tight_layout with padding to accommodate legend
    plt.tight_layout()
    
    # Adjust layout to prevent legend overlap
    plt.subplots_adjust(right=0.85)
    
    plt.savefig(f"{figpath}th_min_map_dbscan_{airstream}_{region}_{timeframe}.png", 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print cluster statistics
    print(f"\nDBSCAN Results:")
    print(f"Number of clusters: {n_clusters}")
    print(f"Number of noise points: {n_noise}")
    print(f"Map extent: {extent}")
    for k in unique_labels:
        if k != -1:
            cluster_data = df_clean[df_clean['cluster'] == k]
            print(f"Cluster {k+1}: {len(cluster_data)} points")
            print(f"  Mean location: ({cluster_data['TH_min_lon'].mean():.2f}°E, {cluster_data['TH_min_lat'].mean():.2f}°N)")
            print(f"  Mean TH_min: {cluster_data['TH_min'].mean():.2f} K")

# %%
def plot_extreme_cooling_scatter_detailed(region, airstream, percentile_threshold=1):
    """
    Detailed scatter plot with regression lines and statistics
    """
    from scipy import stats
    
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 16,
        'legend.fontsize': 12,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14
    })
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    for i, timeframe in enumerate(['present', 'future']):
        ax = axes[i]
        
        # Load and process data
        path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
        df = pd.read_csv(path)
        df_air = df[df['airstream'] == airstream].copy()
        df_clean = df_air.dropna(subset=['cooling_rate_6h', 'cooling_time'])
        
        # Get extreme cooling trajectories
        cooling_threshold = np.percentile(df_clean['cooling_rate_6h'], percentile_threshold)
        extreme_cooling_trajs = df_clean[df_clean['cooling_rate_6h'] <= cooling_threshold]
        
        # Scatter plot
        color = colors[timeframe]
        ax.scatter(extreme_cooling_trajs['cooling_time'], extreme_cooling_trajs['cooling_rate_6h'], 
                  c=color, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
        
        # Add regression line
        if len(extreme_cooling_trajs) > 2:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                extreme_cooling_trajs['cooling_time'], extreme_cooling_trajs['cooling_rate_6h'])
            
            x_reg = np.linspace(extreme_cooling_trajs['cooling_time'].min(), 
                              extreme_cooling_trajs['cooling_time'].max(), 100)
            y_reg = slope * x_reg + intercept
            ax.plot(x_reg, y_reg, '--', color='red', linewidth=2, alpha=0.8)
            
            # Add statistics text
            stats_text = f'R² = {r_value**2:.3f}\np = {p_value:.3f}\nN = {len(extreme_cooling_trajs)}\nThreshold = {cooling_threshold:.3f}'
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Customize subplot
        ax.set_xlabel('Cooling Time [h]')
        ax.set_ylabel('Cooling Rate [K 6h$^{-1}$]')
        ax.set_title(f'{timeframe.capitalize()} - {airstream} - {region}')
        # Set x-axis ticks every 12 hours
        min_time = extreme_cooling_trajs['cooling_time'].min()
        max_time = extreme_cooling_trajs['cooling_time'].max()
        x_min = np.floor(min_time / 12) * 12
        x_max = np.ceil(max_time / 12) * 12
        ax.set_xticks(np.arange(x_min, x_max + 12, 12))
        ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Extreme Cooling Trajectories (Bottom {percentile_threshold}%)', fontsize=18, y=1.02)
    plt.tight_layout()
    
    # Save
    plt.savefig(f"{figpath}extreme_cooling_scatter_detailed_{airstream}_{region}.png", 
                dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}extreme_cooling_scatter_detailed_{airstream}_{region}.pdf", 
                bbox_inches='tight')
    plt.show()
    
    return fig, axes
def plot_cooling_time_histogram_filtered(region, airstream, percentile_threshold=1, bin_width=6):
    """
    Create histogram of cooling times for extreme cooling trajectories
    using the same filtering as the original function
    """
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 16,
        'legend.fontsize': 12,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14
    })
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    for i, timeframe in enumerate(['present', 'future']):
        ax = axes[i]
        
        # Load and process data (same as original function)
        path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
        df = pd.read_csv(path)
        df_air = df[df['airstream'] == airstream].copy()
        df_clean = df_air.dropna(subset=['cooling_rate_6h', 'cooling_time'])
        
        # Get extreme cooling trajectories (same threshold calculation)
        cooling_threshold = np.nanpercentile(df_clean['cooling_rate_6h'].to_numpy(), percentile_threshold)
        extreme_cooling_trajs = df_clean[df_clean['cooling_rate_6h'] <= cooling_threshold]
        
        # Extract cooling times
        cooling_times = extreme_cooling_trajs['cooling_time']
        
        if len(cooling_times) == 0:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
            ax.set_title(f'{timeframe.capitalize()} - No Data')
            continue
        
        # Create bins with 6-hour width
        min_time = cooling_times.min()
        max_time = cooling_times.max()
        bins = np.arange(min_time, max_time + bin_width, bin_width)
        
        # Create histogram
        color = colors[timeframe] if 'colors' in globals() else ('blue' if timeframe == 'present' else 'red')
        n, bin_edges, patches = ax.hist(cooling_times, bins=bins, alpha=0.7, 
                                       color=color, edgecolor='black', linewidth=1)
        
        # Calculate statistics
        mean_time = cooling_times.mean()
        median_time = cooling_times.median()
        std_time = cooling_times.std()
        
        # Add vertical lines for mean and median
        ax.axvline(mean_time, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {mean_time:.1f}h')
        ax.axvline(median_time, color='orange', linestyle='--', linewidth=2, 
                   label=f'Median: {median_time:.1f}h')
        
        # Add statistics text box
        stats_text = f'N = {len(cooling_times)}\nMean = {mean_time:.1f}h\nMedian = {median_time:.1f}h\nStd = {std_time:.1f}h\nThreshold = {cooling_threshold:.3f}K'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, 
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
               fontsize=10)
        
        # Customize subplot with 12-hourly x-axis ticks
        ax.set_xlabel('Cooling Time [h]')
        ax.set_ylabel('Frequency')
        ax.set_title(f'{timeframe.capitalize()} - {airstream} - {region}')
        
        # Set x-axis ticks every 12 hours
        x_min = np.floor(min_time / 12) * 12
        x_max = np.ceil(max_time / 12) * 12
        ax.set_xticks(np.arange(x_min, x_max + 12, 12))
        
        ax.grid(True, alpha=0.3)
        ax.legend(loc='lower right')
        
        # Print statistics for this timeframe
        print(f"\n{timeframe.capitalize()} - {airstream} - {region} Cooling Time Statistics:")
        print(f"N = {len(cooling_times)}")
        print(f"Min: {cooling_times.min():.1f}h")
        print(f"Max: {cooling_times.max():.1f}h")
        print(f"Mean: {mean_time:.1f}h")
        print(f"Median: {median_time:.1f}h")
        print(f"25th percentile: {cooling_times.quantile(0.25):.1f}h")
        print(f"75th percentile: {cooling_times.quantile(0.75):.1f}h")
        print(f"Cooling threshold: {cooling_threshold:.3f}K")
    
    plt.suptitle(f'Cooling Time Distribution - Extreme Cooling Trajectories (Bottom {percentile_threshold}%)', 
                 fontsize=18, y=1.02)
    plt.tight_layout()
    
    # Save plots
    plt.savefig(f"{figpath}cooling_time_histogram_{airstream}_{region}.png", 
                dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}cooling_time_histogram_{airstream}_{region}.pdf", 
                bbox_inches='tight')
    plt.show()
    
    return fig, axes

def plot_combined_evolution_rates_cooling_warming_min24h(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    
    Filters out cooling times shorter than 24 hours.
    """

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4      # Y-axis tick font size
})
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    
    # Filter out cooling times shorter than 24 hours
    df_air_present = df_air_present[df_air_present['cooling_time'] >= 24]
    df_air_future = df_air_future[df_air_future['cooling_time'] >= 24]
    
    # Set up the figure with custom width ratios
    pt = 1./72.27
    textwidth = 418.25368
    fig_width = textwidth * pt  # Use textwidth for overall figure width
    fourto3 = 4/3
    figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
    my_width = 6.69291
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    fig, axes = plt.subplots(1, 3, figsize=(my_width, 2.2), 
                            gridspec_kw={'width_ratios': [3, 1, 1,]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First subplot: Evolution plot
    ax1 = axes[0]
    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})', color=colors[timeframe], linewidth=0.8)
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
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--', alpha=0.7, linewidth=0.8)
    
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

    is_sig_cooling, p_value_cooling, t_stat_cooling = compute_significance_ttest(cooling_rate_present.values, cooling_rate_future.values)
    data_cooling = [cooling_rate_present, cooling_rate_future]

    asterisk_cooling = '*' if is_sig_cooling else ''

    # Vertical boxplots
    box1 = ax2.boxplot([data_cooling[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_cooling[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels([f'present\n mean = {data_cooling[0].mean():.2f}\n median = {data_cooling[0].median():.2f}', f'future\n mean{asterisk_cooling} = {data_cooling[1].mean():.2f}\n median = {data_cooling[1].median():.2f}'], fontsize=6)
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)

    # Fifth subplot: Warming rates boxplot
    ax5 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels([f'present\n mean = {data_warming[0].mean():.2f}\n median = {data_warming[0].median():.2f}', f'future\n mean{asterisk_warming} = {data_warming[1].mean():.2f}\n median = {data_warming[1].median():.2f}'], fontsize=6)
    ax5.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_min24h.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_min24h.pdf", bbox_inches='tight')
    plt.show()

def plot_combined_evolution_rates_cooling_warming_no_outliers(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    
    Filters out the top and bottom 1% of cooling rates.
    """

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4      # Y-axis tick font size
})
    # Load data
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    
    # Load trajectory data for evolution plot
    dfs = {}
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
    
    # Filter out top and bottom 1% of cooling rates
    if varname == 'TH':
        cooling_col = 'cooling_rate_6h'
    elif varname == 'Q':
        cooling_col = 'q_rate_6h_cooling'
    
    p1_present = df_air_present[cooling_col].quantile(0.01)
    p99_present = df_air_present[cooling_col].quantile(0.99)
    p1_future = df_air_future[cooling_col].quantile(0.01)
    p99_future = df_air_future[cooling_col].quantile(0.99)
    
    df_air_present = df_air_present[(df_air_present[cooling_col] >= p1_present) & (df_air_present[cooling_col] <= p99_present)]
    df_air_future = df_air_future[(df_air_future[cooling_col] >= p1_future) & (df_air_future[cooling_col] <= p99_future)]
    
    # Set up the figure with custom width ratios
    pt = 1./72.27
    textwidth = 418.25368
    fig_width = textwidth * pt  # Use textwidth for overall figure width
    fourto3 = 4/3
    figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}
    my_width = 6.69291
    
    # Create subplots with width ratios: evolution plot gets 3 units, each boxplot gets 1 unit
    fig, axes = plt.subplots(1, 3, figsize=(my_width, 2.2), 
                            gridspec_kw={'width_ratios': [3, 1, 1,]})
    
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")
    
    # First subplot: Evolution plot
    ax1 = axes[0]
    for timeframe, (rel_times, means, q25s, q75s) in stats_dict.items():
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})', color=colors[timeframe], linewidth=0.8)
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
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--', alpha=0.7, linewidth=0.8)
    
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

    is_sig_cooling, p_value_cooling, t_stat_cooling = compute_significance_ttest(cooling_rate_present.values, cooling_rate_future.values)
    data_cooling = [cooling_rate_present, cooling_rate_future]

    asterisk_cooling = '*' if is_sig_cooling else ''

    # Vertical boxplots
    box1 = ax2.boxplot([data_cooling[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_cooling[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    ax2.set_xticks([1, 2])
    ax2.set_xticklabels([f'present\n mean = {data_cooling[0].mean():.2f}\n median = {data_cooling[0].median():.2f}', f'future\n mean{asterisk_cooling} = {data_cooling[1].mean():.2f}\n median = {data_cooling[1].median():.2f}'], fontsize=6)
    ax2.set_ylabel(rate_label)
    ax2.grid(True, axis='y', alpha=0.3)

    # Fifth subplot: Warming rates boxplot
    ax5 = axes[2]
    if varname == 'TH':
        warming_rate_present = df_air_present['warming_rate_6h'].dropna()
        warming_rate_future = df_air_future['warming_rate_6h'].dropna()
    elif varname == 'Q':
        warming_rate_present = df_air_present['q_rate_6h_warming'].dropna()
        warming_rate_future = df_air_future['q_rate_6h_warming'].dropna()
    data_warming = [warming_rate_present, warming_rate_future]
    is_sig_warming, p_value_warming, t_stat_warming = compute_significance_ttest(warming_rate_present.values, warming_rate_future.values)
    asterisk_warming = '*' if is_sig_warming else ''
    box5 = ax5.boxplot([data_warming[0]], positions=[1], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
                    medianprops=dict(color='blue', linewidth=0.8),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0.8),
                    whiskerprops=dict(color='orange'),
                    capprops=dict(color='orange'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    ax5.set_xticks([1, 2])
    ax5.set_xticklabels([f'present\n mean = {data_warming[0].mean():.2f}\n median = {data_warming[0].median():.2f}', f'future\n mean{asterisk_warming} = {data_warming[1].mean():.2f}\n median = {data_warming[1].median():.2f}'], fontsize=6)
    ax5.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_no_outliers.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}combined_evolution_rates_{varname}_{airstream}_{region}_coolingwarming_no_outliers.pdf", bbox_inches='tight')
    plt.show()

def compare_cooling_rates_fixed_vs_ensemble(region, airstream, timeframe, time_step_start=-144, time_step_end=-138):
    """
    Compare cooling rates calculated two ways:
    1. Direct: cooling rate between two specific time steps across all trajectories
    2. Ensemble: rate derived from ensemble means at those same time steps
    
    Parameters:
    -----------
    region : str
        'JS' or 'NA'
    airstream : str
        Airstream identifier (e.g., 'LL', 'UL')
    timeframe : str
        'present' or 'future'
    time_step_start : float
        Starting relative time in hours (default -144)
    time_step_end : float
        Ending relative time in hours (default -138)
    """
    
    # Load trajectory data
    path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
    df = pd.read_csv(path)
    
    # Filter for specific airstream
    df_air = df[df['airstream'] == airstream].copy()
    
    if df_air.empty:
        print(f"No data for {airstream} in {region} ({timeframe})")
        return
    
    # Method 1: Direct calculation from individual trajectories
    direct_rates = []
    
    for traj_id, group in df_air.groupby('traj_id'):
        group = group.sort_values('rel_time')
        
        # Get TH values at both time steps
        th_start = group.loc[np.isclose(group['rel_time'], time_step_start), 'TH']
        th_end = group.loc[np.isclose(group['rel_time'], time_step_end), 'TH']
        
        if not th_start.empty and not th_end.empty:
            delta_th = th_end.values[0] - th_start.values[0]
            delta_t = time_step_end - time_step_start  # hours
            cooling_rate = delta_th   # K per 6 hour
            direct_rates.append(cooling_rate)
    
    direct_rates = np.array(direct_rates)
    
    # Method 2: Calculate from ensemble means
    rel_times, means, q25s, q75s = calculate_weighted_stats(df_air, 'TH')
    
    # Find ensemble values at the two time steps (interpolate if necessary)
    th_mean_start = np.interp(time_step_start, rel_times, means)
    th_mean_end = np.interp(time_step_end, rel_times, means)
    delta_th_ensemble = th_mean_end - th_mean_start
    delta_t = time_step_end - time_step_start
    ensemble_rate = delta_th_ensemble   # K per 6 hour
    
    # Statistics
    print(f"\n{'='*70}")
    print(f"Cooling Rate Comparison: {region} - {airstream} - {timeframe}")
    print(f"Time interval: {time_step_start}h to {time_step_end}h (Δt = {delta_t}h)")
    print(f"{'='*70}")
    
    print(f"\nDirect Method (individual trajectories):")
    print(f"  Number of trajectories with data: {len(direct_rates)}")
    print(f"  Mean cooling rate: {direct_rates.mean():.6f} K/6h")
    print(f"  Median cooling rate: {np.median(direct_rates):.6f} K/6h")
    print(f"  Std dev: {direct_rates.std():.6f} K/6h")
    print(f"  Min: {direct_rates.min():.6f} K/6h")
    print(f"  Max: {direct_rates.max():.6f} K/6h")
    print(f"  25th percentile: {np.percentile(direct_rates, 25):.6f} K/6h")
    print(f"  75th percentile: {np.percentile(direct_rates, 75):.6f} K/6h")
    
    print(f"\nEnsemble Method (from means):")
    print(f"  TH mean at {time_step_start}h: {th_mean_start:.4f} K")
    print(f"  TH mean at {time_step_end}h: {th_mean_end:.4f} K")
    print(f"  Ensemble cooling rate: {ensemble_rate:.6f} K/6h")
    
    print(f"\nDifference:")
    print(f"  Direct mean - Ensemble: {direct_rates.mean() - ensemble_rate:.6f} K/6h")
    print(f"  Relative difference: {(direct_rates.mean() - ensemble_rate)/abs(ensemble_rate)*100:.2f}%")
    
    # Create comparison visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram of direct rates
    axes[0].hist(direct_rates, bins=30, alpha=0.7, color='blue', edgecolor='black')
    axes[0].axvline(direct_rates.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {direct_rates.mean():.4f} K/6h')
    axes[0].axvline(ensemble_rate, color='orange', linestyle='--', linewidth=2, label=f'Ensemble: {ensemble_rate:.4f} K/6h')
    axes[0].set_xlabel('Cooling Rate [K/6h]')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title(f'Distribution of Direct Cooling Rates\n{region} - {airstream} - {timeframe}')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Ensemble time series with highlighted time step
    axes[1].plot(rel_times, means, 'b-', linewidth=2, label='Ensemble mean')
    axes[1].fill_between(rel_times, q25s, q75s, alpha=0.2, color='blue', label='IQR')
    
    # Mark the two time steps
    axes[1].plot([time_step_start, time_step_end], [th_mean_start, th_mean_end], 'ro-', 
                linewidth=2, markersize=8, label=f'Rate interval')
    
    # Add slope line
    axes[1].plot([time_step_start, time_step_end], [th_mean_start, th_mean_end], 'r--', 
                alpha=0.5, linewidth=1)
    
    axes[1].set_xlabel('Relative time [h]')
    axes[1].set_ylabel('TH [K]')
    axes[1].set_title(f'Ensemble Mean Temperature Evolution\n{region} - {airstream} - {timeframe}')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xlim(-168, 0)
    
    plt.tight_layout()
    plt.savefig(f"{figpath}cooling_rate_comparison_{region}_{airstream}_{timeframe}_{time_step_start}to{time_step_end}.png", 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    return {
        'direct_rates': direct_rates,
        'ensemble_rate': ensemble_rate,
        'n_trajectories': len(direct_rates)
    }


def plot_cooling_rate_histogram_percent(region, airstream, timeframe, time_step_start=-96, time_step_end=-90, bin_width=0.5):
    """
    Plot only the direct cooling-rate histogram from ``compare_cooling_rates_fixed_vs_ensemble``
    with bins of ``bin_width`` K/6h and percentages on the y-axis.
    """

    # Load trajectory data
    path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
    df = pd.read_csv(path)

    # Filter for specific airstream
    df_air = df[df['airstream'] == airstream].copy()

    if df_air.empty:
        print(f"No data for {airstream} in {region} ({timeframe})")
        return

    direct_rates = []
    for _, group in df_air.groupby('traj_id'):
        group = group.sort_values('rel_time')
        th_start = group.loc[np.isclose(group['rel_time'], time_step_start), 'TH']
        th_end = group.loc[np.isclose(group['rel_time'], time_step_end), 'TH']
        if not th_start.empty and not th_end.empty:
            delta_th = th_end.values[0] - th_start.values[0]
            direct_rates.append(delta_th)

    direct_rates = np.array(direct_rates)
    if direct_rates.size == 0:
        print(f"No matching time steps for {airstream} in {region} ({timeframe})")
        return

    # Ensemble rate for reference line
    rel_times, means, _, _ = calculate_weighted_stats(df_air, 'TH')
    th_mean_start = np.interp(time_step_start, rel_times, means)
    th_mean_end = np.interp(time_step_end, rel_times, means)
    ensemble_rate = th_mean_end - th_mean_start

    # Build bins with requested width
    min_rate = direct_rates.min()
    min_rate = -11.810000000000002
    print(f"Min direct rate: {min_rate}")
    max_rate = direct_rates.max()
    max_rate = 12.25
    print(f"Max direct rate: {max_rate}")
    start_bin = bin_width * np.floor(min_rate / bin_width)
    end_bin = bin_width * np.ceil(max_rate / bin_width)
    bins = np.arange(start_bin, end_bin + bin_width, bin_width)

    # Convert counts to percentages
    weights = np.ones_like(direct_rates) * (100.0 / len(direct_rates))

    plt.rcParams.update({
        'font.size': 6,           # Default font size for ticks
        'axes.titlesize': 6,      # Title font size
        'axes.labelsize': 6,      # X and Y label font size
        'legend.fontsize': 6,     # Legend font size
        'xtick.labelsize': 4,     # X-axis tick font size
        'ytick.labelsize': 4      # Y-axis tick font size
    })


    my_width = 6.69291
    
    fig, ax = plt.subplots(figsize=(my_width/2, 2.2))
   
    # Print the total figure size
    print(f"Total figure size: {fig.get_size_inches()[0]:.2f} x {fig.get_size_inches()[1]:.2f} inches")
    print(f"Total figure size: {fig.get_size_inches()[0]*2.54:.1f} x {fig.get_size_inches()[1]*2.54:.1f} cm")

    if timeframe == 'present':
        fig_label = 'a'
    elif timeframe == 'future':
        fig_label = 'b'

    ax.hist(direct_rates, bins=bins, weights=weights, alpha=0.7, color='blue', edgecolor='black'
            # , label=f'(n={len(direct_rates)})'
            )
    ax.set_ylim(0, 50)
    ax.axvline(direct_rates.mean(), color='red', linestyle='--', linewidth=0.8, label=f'mean: {direct_rates.mean():.4f} K/6h')
    # ax.axvline(ensemble_rate, color='orange', linestyle='--', linewidth=2, label=f'Ensemble: {ensemble_rate:.4f} K/6h')
    ax.set_xlabel('$\\Theta$ rate [K 6h$^{-1}$]')
    ax.set_ylabel('frequency [%]')
    ax.set_title(
        f'{fig_label}) {region} - {airstream} - {timeframe}'
        # f'\n({time_step_start}h to {time_step_end}h)'
    , loc='left')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{figpath}cooling_rate_hist_percent_{region}_{airstream}_{timeframe}_{time_step_start}to{time_step_end}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}cooling_rate_hist_percent_{region}_{airstream}_{timeframe}_{time_step_start}to{time_step_end}.pdf", bbox_inches='tight')
    plt.show()

    return {
        'direct_rates': direct_rates,
        'ensemble_rate': ensemble_rate,
        'bins': bins
    }

def plot_cooling_time_histogram(region, airstream):
    """
    Plot histogram of cooling times for all trajectories in a region/airstream
    with vertical lines showing the ensemble minimum TH time
    """
    
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 14,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Load trajectory data to calculate ensemble minimum time
    dfs = {}
    for timeframe in ['present', 'future']:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
        dfs[timeframe] = pd.read_csv(path)
    
    for i, timeframe in enumerate(['present', 'future']):
        ax = axes[i]
        
        # Load summary data
        path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
        df = pd.read_csv(path)
        
        # Filter for specific airstream
        df_air = df[df['airstream'] == airstream].copy()
        
        if df_air.empty:
            print(f"No data for {airstream} in {region} ({timeframe})")
            continue
        
        # Calculate ensemble minimum time from trajectories
        df_traj = dfs[timeframe]
        df_traj_air = df_traj[df_traj['airstream'] == airstream]
        rel_times, means, q25s, q75s = calculate_weighted_stats(df_traj_air, 'TH')
        min_th_idx = np.argmin(means)
        ensemble_min_time = rel_times[min_th_idx]
        ensemble_min_th = means[min_th_idx]
        
        # Convert ensemble min time to cooling time (hours from start to minimum)
        ensemble_cooling_time = ensemble_min_time + 168
        
        # Extract cooling times
        cooling_times = df_air['cooling_time'].dropna()
        
        if len(cooling_times) == 0:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
            ax.set_title(f'{timeframe.capitalize()} - No Data')
            continue
        
        # Create histogram with 6-hour bins
        bins = np.arange(cooling_times.min(), cooling_times.max() + 6, 6)
        color = colors[timeframe]
        n, bin_edges, patches = ax.hist(cooling_times, bins=bins, alpha=0.7, 
                                       color=color, edgecolor='black', linewidth=1)
        
        # Calculate statistics
        mean_time = cooling_times.mean()
        median_time = cooling_times.median()
        std_time = cooling_times.std()
        q25 = cooling_times.quantile(0.25)
        q75 = cooling_times.quantile(0.75)
        
        # Add vertical lines
        ax.axvline(mean_time, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {mean_time:.1f}h')
        ax.axvline(median_time, color='darkred', linestyle='-.', linewidth=2, 
                   label=f'Median: {median_time:.1f}h')
        ax.axvline(ensemble_cooling_time, color='green', linestyle=':', linewidth=2.5,
                   label=f'Ensemble min: {ensemble_cooling_time:.1f}h')
        
        # Add statistics text box
        stats_text = (f'N = {len(cooling_times)}\n'
                     f'Mean = {mean_time:.1f}h\n'
                     f'Median = {median_time:.1f}h\n'
                     f'Std = {std_time:.1f}h\n'
                     f'Q25 = {q25:.1f}h\n'
                     f'Q75 = {q75:.1f}h\n'
                     f'Ensemble TH_min: {ensemble_min_th:.2f}K')
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, 
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
               fontsize=10)
        
        # Customize subplot
        ax.set_xlabel('Cooling Time [h]')
        ax.set_ylabel('Frequency')
        ax.set_title(f'{timeframe.capitalize()} - {airstream} - {region}')
        ax.grid(True, alpha=0.3, axis='y')
        ax.legend(loc='upper left')
        
        # Print statistics for this timeframe
        print(f"\n{timeframe.capitalize()} - {airstream} - {region} Cooling Time Statistics:")
        print(f"N = {len(cooling_times)}")
        print(f"Min: {cooling_times.min():.1f}h")
        print(f"Max: {cooling_times.max():.1f}h")
        print(f"Mean: {mean_time:.1f}h")
        print(f"Median: {median_time:.1f}h")
        print(f"Std: {std_time:.1f}h")
        print(f"Q25: {q25:.1f}h")
        print(f"Q75: {q75:.1f}h")
        print(f"Ensemble minimum TH time: {ensemble_min_time:.1f}h (cooling time: {ensemble_cooling_time:.1f}h)")
        print(f"Ensemble minimum TH value: {ensemble_min_th:.2f}K")
    
    plt.suptitle(f'Cooling Time Distribution: {airstream} - {region}', fontsize=16, y=1.00)
    plt.tight_layout()
    
    # Save plots
    plt.savefig(f"{figpath}cooling_time_histogram_{airstream}_{region}.png", 
                dpi=300, bbox_inches='tight')
    plt.savefig(f"{figpath}cooling_time_histogram_{airstream}_{region}.pdf", 
                bbox_inches='tight')
    plt.show()
    
    return fig

def count_cooling_time(region, airstream=None, cooling_time_value=0):
    """
    Count trajectories with cooling_time == 0 for a given region and optionally airstream
    
    Parameters:
    -----------
    region : str
        'JS' or 'NA'
    airstream : str, optional
        Airstream identifier (e.g., 'LL', 'UL'). If None, checks all airstreams.
    """
    
    for timeframe in ['present', 'future']:
        path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
        df = pd.read_csv(path)
        
        if airstream is not None:
            df = df[df['airstream'] == airstream]
        
        # Count trajectories with cooling_time == cooling_time_value
        zero_cooling = df[df['cooling_time'] == cooling_time_value]
        total = len(df)
        
        print(f"\n{region} - {timeframe}" + (f" - {airstream}" if airstream else " - ALL"))
        print(f"  Trajectories with cooling_time = {cooling_time_value}: {len(zero_cooling)} / {total} ({len(zero_cooling)/total*100:.2f}%)")
        
        if len(zero_cooling) > 0:
            print(f"  Their TH_min_time values: {zero_cooling['TH_min_time'].unique()}")
            print(f"  Their TH_min values range: [{zero_cooling['TH_min'].min():.2f}, {zero_cooling['TH_min'].max():.2f}] K")


def cooling_rate_distribution(region, airstream, timestep_start=-144, timestep_end=-138):
    """
    Plot distribution of cooling rates for given region and airstream and timestep.
    
    :param region: Description
    :param airstream: Description
    """
    timestep_end = timestep_start + 6  # Ensure 6-hour interval
    rates_dict = {}

    for timeframe in ['present', 'future']:
        # Load trajectory data
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
        df = pd.read_csv(path)
        
        # Filter for specific airstream
        df_air = df[df['airstream'] == airstream].copy()
        
        if df_air.empty:
            print(f"No data for {airstream} in {region} ({timeframe})")
            return
        
        # Method 1: Direct calculation from individual trajectories
        direct_rates = []
        
        for traj_id, group in df_air.groupby('traj_id'):
            group = group.sort_values('rel_time')
            
            # Get TH values at both time steps
            th_start = group.loc[np.isclose(group['rel_time'], timestep_start), 'TH']
            th_end = group.loc[np.isclose(group['rel_time'], timestep_end), 'TH']
            
            if not th_start.empty and not th_end.empty:
                delta_th = th_end.values[0] - th_start.values[0]
                delta_t = timestep_end - timestep_start  # hours
                cooling_rate = delta_th   # K per 6 hours
                direct_rates.append(cooling_rate)

        # count negative and positive rates
        direct_rates = np.array(direct_rates)
        n_negative = np.sum(direct_rates < 0)
        n_positive = np.sum(direct_rates > 0)
        n_zero = np.sum(direct_rates == 0)
        print(f"\n{region} - {airstream} - {timeframe} Cooling Rate Distribution ({timestep_start}h to {timestep_end}h):")
        print(f"  Total trajectories: {len(direct_rates)}")
        print(f"  Negative rates: {n_negative} ({n_negative/len(direct_rates)*100:.2f}%)")
        print(f"  Positive rates: {n_positive} ({n_positive/len(direct_rates)*100:.2f}%)")
        print(f"  Zero rates: {n_zero} ({n_zero/len(direct_rates)*100:.2f}%)")

        # calculate mean of negative and positive rates separately
        mean_negative = direct_rates[direct_rates < 0].mean() if n_negative > 0 else np.nan
        mean_positive = direct_rates[direct_rates > 0].mean() if n_positive > 0 else np.nan
        print(f"  Mean negative rate: {mean_negative:.4f} K/6h")
        print(f"  Mean positive rate: {mean_positive:.4f} K/6h")

        rates_dict[timeframe] = np.array(direct_rates)
    
    # Plotting
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 16,
        'axes.labelsize': 14,
        'legend.fontsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })
    fig, ax = plt.subplots(figsize=(8, 6))
    bins = np.linspace(-12, 12, 100)
    for timeframe, rates in rates_dict.items():
        ax.hist(rates, bins=bins, alpha=0.5, label=f'{timeframe.capitalize()} (N={len(rates)})', edgecolor='black')
        ax.axvline(x=rates.mean(), color='blue' if timeframe=='present' else 'orange', linestyle='--', 
                  label=f'Mean {timeframe}: {rates.mean():.2f} K/6h')
    ax.set_xlabel('Cooling Rate [K/6h]')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Cooling Rate Distribution: {airstream} - {region}\n({timestep_start}h to {timestep_end}h)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{figpath}cooling_rate_distribution_{region}_{airstream}_{timestep_start}to{timestep_end}.png", dpi=300, bbox_inches='tight')
    plt.show()
# # %%

# # --- First pass: determine auto y-limits for each REGION ---
# auto_limits = {region: [] for region in regions}  # Changed from airstream to region
# varname = 'TH'  # Change this to analyze a different variable
# rounding = 5  # round y-limits to the nearest digit (5 for TH, 1 for Q, 50 for HFLX)

# for region in regions:
#     dfs = {}
#     for timeframe in timeframes:
#         path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
#         dfs[timeframe] = pd.read_csv(path)
#     for airstream in airstreams:
#         fig, ax = plt.subplots()
#         for timeframe in timeframes:
#             df = dfs[timeframe]
#             df_air = df[df['airstream'] == airstream]
#             if df_air.empty:
#                 continue
#             rel_times, means, q25s, q75s = calculate_weighted_stats(df_air, varname)
#             ax.plot(rel_times, means, color=colors[timeframe])
#             ax.fill_between(rel_times, q25s, q75s, alpha=0.2, color=colors[timeframe])
#         # Record the auto y-limits for the REGION
#         auto_limits[region].append(ax.get_ylim())  # Changed from airstream to region
#         plt.close(fig)

# # Find the global min/max of the auto-limits for each REGION, round to X
# th_limits = {}
# for region in regions:  # Changed from airstream to region
#     if auto_limits[region]:  # Changed from airstream to region
#         ymins, ymaxs = zip(*auto_limits[region])  # Changed from airstream to region
#         th_min = rounding * np.floor(min(ymins) / rounding)
#         th_max = rounding * np.ceil(max(ymaxs) / rounding)
#         th_limits[region] = (th_min, th_max)  # Changed from airstream to region
#     else:
#         # Default limits if no data available
#         th_limits[region] = (250, 350)  # Changed from airstream to region
# if varname == 'Q':
#     for region in regions:
#         th_limits[region] = (0, 5)
# # --- Second pass: plot with consistent y-limits ---
# for region in regions:
#     dfs = {}
#     for timeframe in timeframes:
#         path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
#         dfs[timeframe] = pd.read_csv(path)
#     for airstream in airstreams:
#         stats_dict = {}
#         for timeframe in timeframes:
#             df = dfs[timeframe]
#             df_air = df[df['airstream'] == airstream]
#             if df_air.empty:
#                 continue
#             stats_dict[timeframe] = calculate_weighted_stats(df_air, varname)

#             if varname == 'TH':
#                 min_idx = np.argmin(stats_dict[timeframe][1])
#                 ensemble_min_time = stats_dict[timeframe][0][min_idx]
#                 ensemble_min_th = stats_dict[timeframe][1][min_idx]
#                 ensemble_168_th = stats_dict[timeframe][1][0]
#                 ensemble_0_th = stats_dict[timeframe][1][-1]
                
#                 cooling = ensemble_min_th - ensemble_168_th
#                 warming = ensemble_0_th - ensemble_min_th

#                 cooling_time = ensemble_min_time + 168
#                 warming_time = 0 - ensemble_min_time

#                 cooling_slope = cooling/cooling_time
#                 warming_slope = warming/warming_time

#                 print(f"Stats for {airstream} in {region} ({timeframe}):")
#                 print(f'Number of trajectories: {len(df_air["traj_id"].unique())}')
#                 print(f'Ensemble min time: {ensemble_min_time} h')
#                 print(f'Ensemble min TH: {ensemble_min_th:.2f} K')
#                 print(f'TH at -168 h: {ensemble_168_th:.2f} K')
#                 print(f'TH at 0 h: {ensemble_0_th:.2f} K')
#                 print(f'Cooling: {cooling:.2f} K, Warming: {warming:.2f} K')
#                 print(f'Cooling time: {cooling_time:.2f} h, Warming time: {warming_time:.2f} h')
#                 print(f'Cooling slope: {cooling_slope:.2f} K/h, Warming slope: {warming_slope:.2f} K/h')
#                 print(f'Cooling slope: {cooling_slope*24:.2f} K/d, Warming slope: {warming_slope*24:.2f} K/d')
#                 print('--------------------------------------------')
#         if stats_dict:
#             plot_weighted_stats(stats_dict, th_limits, region, airstream, varname)
# %% NEW
# For each traj_id, find TH at rel_time=0, TH at rel_time=-168, and minimum TH
regions = ['JS', 'NA']
timeframes = ['present', 'future']
airstreams = ['LL', 'UL']

# regions = ['NA']  # For testing, can be expanded later
# timeframes = ['present']
# airstreams = ['LL']

for region in regions:
    for timeframe in timeframes:
        path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
                land_mask = (group['rel_time'] > TH_min_time) & (group['rel_time'] <= 0) & ~group['TH_6h_change_rate'].isna() & (group['SST'] <= 0)

                warming_rates = group.loc[warming_mask, 'TH_6h_change_rate']
                q_rates_warming = group.loc[warming_mask, 'Q_6h_change_rate']
                warming_rate_6h = warming_rates.mean() if len(warming_rates) > 0 else np.nan
                q_rate_6h_warming = q_rates_warming.mean() if len(q_rates_warming) > 0 else np.nan

                warming_rates_land = group.loc[land_mask, 'TH_6h_change_rate']
                warming_rate_6h_land = warming_rates_land.mean() if len(warming_rates_land) > 0 else np.nan
                q_rates_land = group.loc[land_mask, 'Q_6h_change_rate']
                q_rate_6h_land = q_rates_land.mean() if len(q_rates_land) > 0 else np.nan

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
                    'q_rate_6h_wbc': q_rate_6h_wbc,    # NEW: Average Q per 6-hour interval during warming (WBC)
                    'warming_rate_6h_land': warming_rate_6h_land,   # NEW: Average K per 6-hour interval during warming (land)
                    'q_rate_6h_land': q_rate_6h_land,    # NEW: Average Q per 6-hour interval during warming (land)
                })

        summary_df = pd.DataFrame(summary_list)

        # for each traj_id, get the slope of cooling and warming
        summary_df['cooling'] = summary_df['TH_min'] - summary_df['TH_at_-168']
        summary_df['warming'] = summary_df['TH_at_0'] - summary_df['TH_min']
        summary_df['cooling_time'] = summary_df['TH_min_time'] + 168
        summary_df['warming_time'] = 0-summary_df['TH_min_time']

        # Example: save to CSV
        summary_df.to_csv(f"{datapath}TH_summary_{region}_{timeframe}.csv", index=False, header=True)
# %%
# # For each traj_id, find TH at rel_time=0, TH at rel_time=-168, and minimum TH
# regions = ['JS', 'NA']
# timeframes = ['present', 'future']
# airstreams = ['BCAR', 'DI', 'nonBCAR', 'LL']
# airstreams = ['LL', 'UL']

# for region in regions:
#     for timeframe in timeframes:
#         path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
#         df = pd.read_csv(path)
#         summary_list = []
#         for airstream in airstreams:
#             df_air = df[df['airstream'] == airstream]
#             if df_air.empty:
#                 continue
#             grouped = df_air.groupby('traj_id')
#             for traj_id, group in grouped:
#                 th_0 = group.loc[np.isclose(group['rel_time'], 0.0), 'TH']
#                 th_168 = group.loc[np.isclose(group['rel_time'], -168.0), 'TH']
#                 th_min = group['TH'].min()
#                 TH_min_time = group.loc[group['TH'] == th_min, 'rel_time'].values[0] if not group[group['TH'] == th_min].empty else np.nan
#                 summary_list.append({
#                     'region': region,
#                     'timeframe': timeframe,
#                     'airstream': airstream,
#                     'traj_id': traj_id,
#                     'TH_at_0': th_0.iloc[0] if not th_0.empty else np.nan,
#                     'TH_at_-168': th_168.iloc[0] if not th_168.empty else np.nan,
#                     'TH_min': th_min,
#                     'TH_min_time': TH_min_time
#                 })

#         summary_df = pd.DataFrame(summary_list)

#         # for each traj_id, get the slope of cooling and warming
#         summary_df['cooling'] = summary_df['TH_min'] - summary_df['TH_at_-168']
#         summary_df['warming'] = summary_df['TH_at_0'] - summary_df['TH_min']
#         summary_df['cooling_time'] = summary_df['TH_min_time'] + 168
#         summary_df['warming_time'] = 0-summary_df['TH_min_time']
#         summary_df['slope_cooling'] = (summary_df['TH_at_-168'] - summary_df['TH_min']) / (-168 - summary_df['TH_min_time'])
#         summary_df['slope_warming'] = (summary_df['TH_at_0'] - summary_df['TH_min']) / (0 - summary_df['TH_min_time'])
#         summary_df['slope_cooling_d'] = summary_df['slope_cooling'] * 24
#         summary_df['slope_warming_d'] = summary_df['slope_warming'] * 24

#         # Example: save to CSV
#         summary_df.to_csv(f"{datapath}TH_summary_{region}_{timeframe}.csv", index=False, header=True)

# %% NEW 
# Now we have a summary CSV for each region and timeframe with the necessary statistics
# For each region and airstream, make separate boxplots for cooling and warming for present and future

regions = ['JS', 'NA']
airstreams = ['LL', 'UL']

regions = ['NA']  # For testing, can be expanded later
airstreams = ['LL']

plt.rcParams.update({
    'font.size': 16,           # Default font size for ticks
    'axes.titlesize': 20,      # Title font size
    'axes.labelsize': 20,      # X and Y label font size
    'legend.fontsize': 20,     # Legend font size
    'xtick.labelsize': 16,     # X-axis tick font size
    'ytick.labelsize': 16      # Y-axis tick font size
})

for region in regions:
    print('--------------------------------------------')
    print('--------------------------------------------')
    print(f"Region: {region}")
    print('--------------------------------------------')
    path_present = f"{datapath}TH_summary_{region}_present.csv"
    df_present = pd.read_csv(path_present)
    path_future = f"{datapath}TH_summary_{region}_future.csv"
    df_future = pd.read_csv(path_future)
    for airstream in airstreams:
        df_air_present = df_present[df_present['airstream'] == airstream]
        df_air_future = df_future[df_future['airstream'] == airstream]
        if df_air_present.empty or df_air_future.empty:
            continue

#         #mask out zero values and nans
#         cooling_present = df_air_present['cooling'].where(df_air_present['cooling']<0).dropna()
#         warming_present = df_air_present['warming'].where(df_air_present['warming']>0).dropna()
#         cooling_future = df_air_future['cooling'].where(df_air_future['cooling']<0).dropna()
#         warming_future = df_air_future['warming'].where(df_air_future['warming']>0).dropna()

#         cooling_rate_6h_present = df_air_present['cooling_rate_6h'].dropna()
#         warming_rate_6h_present = df_air_present['warming_rate_6h'].dropna()
#         cooling_rate_6h_future = df_air_future['cooling_rate_6h'].dropna()
#         warming_rate_6h_future = df_air_future['warming_rate_6h'].dropna()

#         print('--------------------------------------------')
#         print(f"Airstream: {airstream}")
#         print('data sizes:')
#         print(f'cooling rate present: {len(cooling_rate_6h_present)}, cooling rate future: {len(cooling_rate_6h_future)}')
#         print(f'warming rate present: {len(warming_rate_6h_present)}, warming rate future: {len(warming_rate_6h_future)}')
#         print('Cooling rate significance (present vs. future):')
#         print(compute_significance_ttest(cooling_rate_6h_present, cooling_rate_6h_future))
#         print(mannwhitneyu(cooling_rate_6h_present, cooling_rate_6h_future))
#         print('Warming rate significance (present vs. future):')
#         print(compute_significance_ttest(warming_rate_6h_present, warming_rate_6h_future))
#         print(mannwhitneyu(warming_rate_6h_present, warming_rate_6h_future))

#         # Cooling: histogram + horizontal boxplot
#         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
#         # fig, axes = plt.subplots(2, 1, figsize=(figsizes["twocol"], figsizes["twocol"]), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)

        
#         # Histogram
#         axes[0].hist(cooling_present, bins=50, alpha=0.5, label='present', color='blue')
#         axes[0].hist(cooling_future, bins=50, alpha=0.5, label='future', color='orange')
#         axes[0].set_ylabel('frequency')
#         # axes[0].set_title(f'{airstream} - {region} - Cooling')
#         axes[0].legend()
#         axes[0].grid(True, axis='y', alpha=0.3)
#         # Horizontal boxplot
#         data_cooling = [cooling_present, cooling_future]
#         axes[1].boxplot(data_cooling[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
#                         medianprops=dict(color='blue', linewidth=1.5),
#                         whiskerprops=dict(color='blue'),
#                         capprops=dict(color='blue'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].boxplot([data_cooling[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
#                         medianprops=dict(color='darkorange', linewidth=1.5),
#                         whiskerprops=dict(color='orange'),
#                         capprops=dict(color='orange'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].set_yticks([])
#         # axes[1].set_yticklabels(['present', 'future'])
#         axes[1].set_xlabel('cooling [K]')
#         axes[1].grid(True, axis='x', alpha=0.3)
#         plt.tight_layout()
#         plt.savefig(f"{figpath}hist_box_cooling_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
#         plt.show()

#         # Warming: histogram + horizontal boxplot
#         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
#         axes[0].hist(warming_present, bins=50, alpha=0.5, label='present', color='blue')
#         axes[0].hist(warming_future, bins=50, alpha=0.5, label='future', color='orange')
#         axes[0].set_ylabel('frequency')
#         # axes[0].set_title(f'{airstream} - {region} - Warming')
#         axes[0].legend()
#         axes[0].grid(True, axis='y', alpha=0.3)
#         data_warming = [warming_present, warming_future]
#         axes[1].boxplot(data_warming[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
#                         medianprops=dict(color='blue', linewidth=1.5),
#                         whiskerprops=dict(color='blue'),
#                         capprops=dict(color='blue'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].boxplot([data_warming[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
#                         medianprops=dict(color='darkorange', linewidth=1.5),
#                         whiskerprops=dict(color='orange'),
#                         capprops=dict(color='orange'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].set_yticks([])
#         axes[1].set_xlabel('warming [K]')
#         axes[1].grid(True, axis='x', alpha=0.3)
#         plt.tight_layout()
#         plt.savefig(f"{figpath}hist_box_warming_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
#         plt.show()

#         # Cooling rate (K/6h)
#         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
#         axes[0].hist(cooling_rate_6h_present, bins=50, alpha=0.5, label='present', color='blue')
#         axes[0].hist(cooling_rate_6h_future, bins=50, alpha=0.5, label='future', color='orange')
#         axes[0].set_ylabel('frequency')
#         # axes[0].set_title(f'{airstream} - {region} - Cooling Slope')
#         axes[0].legend()
#         axes[0].grid(True, axis='y', alpha=0.3)
#         data_slopes_cooling = [cooling_rate_6h_present, cooling_rate_6h_future]
#         axes[1].boxplot(data_slopes_cooling[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
#                         medianprops=dict(color='blue', linewidth=1.5),
#                         whiskerprops=dict(color='blue'),
#                         capprops=dict(color='blue'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].boxplot([data_slopes_cooling[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
#                         medianprops=dict(color='darkorange', linewidth=1.5),
#                         whiskerprops=dict(color='orange'),
#                         capprops=dict(color='orange'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].set_yticks([])
#         axes[1].set_xlabel('cooling rate [K 6h$^{-1}$]')
#         if airstream in ['BCAR', 'nonBCAR']:
#             axes[1].set_xlim(-40, 0)
#             axes[1].set_xticks(np.arange(-40, 1, 5))        
#         axes[1].grid(True, axis='x', alpha=0.3)
#         plt.tight_layout()
#         plt.savefig(f"{figpath}hist_box_cooling_slope_6h_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
#         plt.show()

#         # Warming rate (K/6h)
#         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
#         axes[0].hist(warming_rate_6h_present, bins=50, alpha=0.5, label='present', color='blue')
#         axes[0].hist(warming_rate_6h_future, bins=50, alpha=0.5, label='future', color='orange')
#         axes[0].set_ylabel('frequency')
#         # axes[0].set_title(f'{airstream} - {region} - Warming Slope')
#         axes[0].legend()
#         axes[0].grid(True, axis='y', alpha=0.3)
#         data_slopes_warming = [warming_rate_6h_present, warming_rate_6h_future]
#         axes[1].boxplot(data_slopes_warming[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
#                         medianprops=dict(color='blue', linewidth=1.5),
#                         whiskerprops=dict(color='blue'),
#                         capprops=dict(color='blue'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].boxplot([data_slopes_warming[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
#                         medianprops=dict(color='darkorange', linewidth=1.5),
#                         whiskerprops=dict(color='orange'),
#                         capprops=dict(color='orange'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].set_yticks([])
#         axes[1].set_xlabel('warming rate [K 6h$^{-1}$]')
#         if airstream in ['BCAR', 'nonBCAR']:
#             axes[1].set_xlim(0, 80)
#             axes[1].set_xticks(np.arange(0, 81, 10))
#         axes[1].grid(True, axis='x', alpha=0.3)
#         plt.tight_layout()
#         plt.savefig(f"{figpath}hist_box_warming_slope_6h_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
#         plt.show()

#         # Q rate during cooling (Q/6h)
#         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
#         axes[0].hist(df_air_present['q_rate_6h_cooling'].dropna(), bins=50, alpha=0.5, label='present', color='blue')
#         axes[0].hist(df_air_future['q_rate_6h_cooling'].dropna(), bins=50, alpha=0.5, label='future', color='orange')
#         axes[0].set_ylabel('frequency')
#         # axes[0].set_title(f'{airstream} - {region} -
#         axes[0].legend()
#         axes[0].grid(True, axis='y', alpha=0.3)
#         data_q_cooling = [df_air_present['q_rate_6h_cooling'].dropna(), df_air_future['q_rate_6h_cooling'].dropna()]
#         axes[1].boxplot(data_q_cooling[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
#                         medianprops=dict(color='blue', linewidth=1.5),
#                         whiskerprops=dict(color='blue'),
#                         capprops=dict(color='blue'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].boxplot([data_q_cooling[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
#                         medianprops=dict(color='darkorange', linewidth=1.5),
#                         whiskerprops=dict(color='orange'),
#                         capprops=dict(color='orange'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].set_yticks([])
#         axes[1].set_xlabel('Q rate during cooling [g kg$^{-1}$ 6h$^{-1}$]')
#         axes[1].grid(True, axis='x', alpha=0.3)
#         plt.tight_layout()
#         plt.savefig(f"{figpath}hist_box_q_rate_cooling_6h_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
#         plt.show()  

#         # Q rate during warming (Q/6h)
#         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)        
#         axes[0].hist(df_air_present['q_rate_6h_warming'].dropna(), bins=50, alpha=0.5, label='present', color='blue')
#         axes[0].hist(df_air_future['q_rate_6h_warming'].dropna(), bins=50, alpha=0.5, label='future', color='orange')
#         axes[0].set_ylabel('frequency')
#         # axes[0].set_title(f'{airstream} - {region} -
#         axes[0].legend()
#         axes[0].grid(True, axis='y', alpha=0.3)
#         data_q_warming = [df_air_present['q_rate_6h_warming'].dropna(), df_air_future['q_rate_6h_warming'].dropna()]
#         axes[1].boxplot(data_q_warming[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
#                         medianprops=dict(color='blue', linewidth=1.5),
#                         whiskerprops=dict(color='blue'),
#                         capprops=dict(color='blue'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].boxplot([data_q_warming[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
#                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
#                         medianprops=dict(color='darkorange', linewidth=1.5),
#                         whiskerprops=dict(color='orange'),
#                         capprops=dict(color='orange'),
#                         whis=(1, 99),
#                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
#                         showmeans=True,
#                         meanline=True,
#                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
#         axes[1].set_yticks([])
#         axes[1].set_xlabel('Q rate during warming [g kg$^{-1}$ 6h$^{-1}$]')
#         axes[1].grid(True, axis='x', alpha=0.3)
#         plt.tight_layout()
#         plt.savefig(f"{figpath}hist_box_q_rate_warming_6h_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
#         plt.show()

# # # %%
# # # Now we have a summary CSV for each region and timeframe with the necessary statistics
# # # For each region and airstream, make separate boxplots for cooling and warming for present and future

# # regions = ['JS', 'NA']
# # airstreams = ['BCAR', 'DI', 'nonBCAR', 'LL']
# # airstreams = ['LL', 'UL']


# # # regions = ['NA']  # For testing, can be expanded later
# # # airstreams = ['BCAR']

# # plt.rcParams.update({
# #     'font.size': 16,           # Default font size for ticks
# #     'axes.titlesize': 20,      # Title font size
# #     'axes.labelsize': 20,      # X and Y label font size
# #     'legend.fontsize': 20,     # Legend font size
# #     'xtick.labelsize': 16,     # X-axis tick font size
# #     'ytick.labelsize': 16      # Y-axis tick font size
# # })

# # for region in regions:
# #     print('--------------------------------------------')
# #     print('--------------------------------------------')
# #     print(f"Region: {region}")
# #     print('--------------------------------------------')
# #     path_present = f"{datapath}TH_summary_{region}_present.csv"
# #     df_present = pd.read_csv(path_present)
# #     path_future = f"{datapath}TH_summary_{region}_future.csv"
# #     df_future = pd.read_csv(path_future)
# #     for airstream in airstreams:
# #         df_air_present = df_present[df_present['airstream'] == airstream]
# #         df_air_future = df_future[df_future['airstream'] == airstream]
# #         if df_air_present.empty or df_air_future.empty:
# #             continue

# #         #mask out zero values and nans
# #         cooling_present = df_air_present['cooling'].where(df_air_present['cooling']<0).dropna()
# #         warming_present = df_air_present['warming'].where(df_air_present['warming']>0).dropna()
# #         cooling_future = df_air_future['cooling'].where(df_air_future['cooling']<0).dropna()
# #         warming_future = df_air_future['warming'].where(df_air_future['warming']>0).dropna()

# #         cooling_slope_present = df_air_present['slope_cooling'].dropna()
# #         warming_slope_present = df_air_present['slope_warming'].dropna()
# #         cooling_slope_future = df_air_future['slope_cooling'].dropna()
# #         warming_slope_future = df_air_future['slope_warming'].dropna()

# #         cooling_slope_d_present = df_air_present['slope_cooling_d'].dropna()
# #         warming_slope_d_present = df_air_present['slope_warming_d'].dropna()
# #         cooling_slope_d_future = df_air_future['slope_cooling_d'].dropna()
# #         warming_slope_d_future = df_air_future['slope_warming_d'].dropna()

# #         # Cooling: histogram + horizontal boxplot
# #         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
# #         # fig, axes = plt.subplots(2, 1, figsize=(figsizes["twocol"], figsizes["twocol"]), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)

        
# #         # Histogram
# #         axes[0].hist(cooling_present, bins=50, alpha=0.5, label='present', color='blue')
# #         axes[0].hist(cooling_future, bins=50, alpha=0.5, label='future', color='orange')
# #         axes[0].set_ylabel('frequency')
# #         # axes[0].set_title(f'{airstream} - {region} - Cooling')
# #         axes[0].legend()
# #         axes[0].grid(True, axis='y', alpha=0.3)
# #         # Horizontal boxplot
# #         data_cooling = [cooling_present, cooling_future]
# #         axes[1].boxplot(data_cooling[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
# #                         medianprops=dict(color='blue', linewidth=1.5),
# #                         whiskerprops=dict(color='blue'),
# #                         capprops=dict(color='blue'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].boxplot([data_cooling[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
# #                         medianprops=dict(color='darkorange', linewidth=1.5),
# #                         whiskerprops=dict(color='orange'),
# #                         capprops=dict(color='orange'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].set_yticks([])
# #         # axes[1].set_yticklabels(['present', 'future'])
# #         axes[1].set_xlabel('cooling [K]')
# #         axes[1].grid(True, axis='x', alpha=0.3)
# #         plt.tight_layout()
# #         plt.savefig(f"{figpath}hist_box_cooling_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# #         # Warming: histogram + horizontal boxplot
# #         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
# #         axes[0].hist(warming_present, bins=50, alpha=0.5, label='present', color='blue')
# #         axes[0].hist(warming_future, bins=50, alpha=0.5, label='future', color='orange')
# #         axes[0].set_ylabel('frequency')
# #         # axes[0].set_title(f'{airstream} - {region} - Warming')
# #         axes[0].legend()
# #         axes[0].grid(True, axis='y', alpha=0.3)
# #         data_warming = [warming_present, warming_future]
# #         axes[1].boxplot(data_warming[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
# #                         medianprops=dict(color='blue', linewidth=1.5),
# #                         whiskerprops=dict(color='blue'),
# #                         capprops=dict(color='blue'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].boxplot([data_warming[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
# #                         medianprops=dict(color='darkorange', linewidth=1.5),
# #                         whiskerprops=dict(color='orange'),
# #                         capprops=dict(color='orange'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].set_yticks([])
# #         axes[1].set_xlabel('warming [K]')
# #         axes[1].grid(True, axis='x', alpha=0.3)
# #         plt.tight_layout()
# #         plt.savefig(f"{figpath}hist_box_warming_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# #         # Cooling slope: histogram + horizontal boxplot
# #         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
# #         axes[0].hist(cooling_slope_present, bins=50, alpha=0.5, label='present', color='blue')
# #         axes[0].hist(cooling_slope_future, bins=50, alpha=0.5, label='future', color='orange')
# #         axes[0].set_ylabel('frequency')
# #         # axes[0].set_title(f'{airstream} - {region} - Cooling Slope')
# #         axes[0].legend()
# #         axes[0].grid(True, axis='y', alpha=0.3)
# #         data_slopes_cooling = [cooling_slope_present, cooling_slope_future]
# #         axes[1].boxplot(data_slopes_cooling[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
# #                         medianprops=dict(color='blue', linewidth=1.5),
# #                         whiskerprops=dict(color='blue'),
# #                         capprops=dict(color='blue'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].boxplot([data_slopes_cooling[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
# #                         medianprops=dict(color='darkorange', linewidth=1.5),
# #                         whiskerprops=dict(color='orange'),
# #                         capprops=dict(color='orange'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].set_yticks([])
# #         axes[1].set_xlabel('cooling rate [K h$^{-1}$]')
# #         axes[1].grid(True, axis='x', alpha=0.3)
# #         plt.tight_layout()
# #         plt.savefig(f"{figpath}hist_box_cooling_slope_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# #         # Warming slope: histogram + horizontal boxplot
# #         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
# #         axes[0].hist(warming_slope_present, bins=50, alpha=0.5, label='present', color='blue')
# #         axes[0].hist(warming_slope_future, bins=50, alpha=0.5, label='future', color='orange')
# #         axes[0].set_ylabel('frequency')
# #         # axes[0].set_title(f'{airstream} - {region} - Warming Slope')
# #         axes[0].legend()
# #         axes[0].grid(True, axis='y', alpha=0.3)
# #         data_slopes_warming = [warming_slope_present, warming_slope_future]
# #         axes[1].boxplot(data_slopes_warming[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
# #                         medianprops=dict(color='blue', linewidth=1.5),
# #                         whiskerprops=dict(color='blue'),
# #                         capprops=dict(color='blue'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].boxplot([data_slopes_warming[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
# #                         medianprops=dict(color='darkorange', linewidth=1.5),
# #                         whiskerprops=dict(color='orange'),
# #                         capprops=dict(color='orange'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].set_yticks([])
# #         axes[1].set_xlabel('warming rate [K h$^{-1}$]')
# #         axes[1].grid(True, axis='x', alpha=0.3)
# #         plt.tight_layout()
# #         plt.savefig(f"{figpath}hist_box_warming_slope_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# #         # Cooling slope (K/d)
# #         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
# #         axes[0].hist(cooling_slope_d_present, bins=50, alpha=0.5, label='present', color='blue')
# #         axes[0].hist(cooling_slope_d_future, bins=50, alpha=0.5, label='future', color='orange')
# #         axes[0].set_ylabel('frequency')
# #         # axes[0].set_title(f'{airstream} - {region} - Cooling Slope')
# #         axes[0].legend()
# #         axes[0].grid(True, axis='y', alpha=0.3)
# #         data_slopes_cooling = [cooling_slope_d_present, cooling_slope_d_future]
# #         axes[1].boxplot(data_slopes_cooling[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
# #                         medianprops=dict(color='blue', linewidth=1.5),
# #                         whiskerprops=dict(color='blue'),
# #                         capprops=dict(color='blue'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].boxplot([data_slopes_cooling[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
# #                         medianprops=dict(color='darkorange', linewidth=1.5),
# #                         whiskerprops=dict(color='orange'),
# #                         capprops=dict(color='orange'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].set_yticks([])
# #         axes[1].set_xlabel('cooling rate [K d$^{-1}$]')
# #         if airstream in ['BCAR', 'nonBCAR']:
# #             axes[1].set_xlim(-40, 0)
# #             axes[1].set_xticks(np.arange(-40, 1, 5))        
# #         axes[1].grid(True, axis='x', alpha=0.3)
# #         plt.tight_layout()
# #         plt.savefig(f"{figpath}hist_box_cooling_slope_d_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# #         # Warming slope (K/d)
# #         fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 2]}, sharex=True)
# #         axes[0].hist(warming_slope_d_present, bins=50, alpha=0.5, label='present', color='blue')
# #         axes[0].hist(warming_slope_d_future, bins=50, alpha=0.5, label='future', color='orange')
# #         axes[0].set_ylabel('frequency')
# #         # axes[0].set_title(f'{airstream} - {region} - Warming Slope')
# #         axes[0].legend()
# #         axes[0].grid(True, axis='y', alpha=0.3)
# #         data_slopes_warming = [warming_slope_d_present, warming_slope_d_future]
# #         axes[1].boxplot(data_slopes_warming[0], vert=False, positions=[1], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='blue', color='blue', alpha=0.5),
# #                         medianprops=dict(color='blue', linewidth=1.5),
# #                         whiskerprops=dict(color='blue'),
# #                         capprops=dict(color='blue'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].boxplot([data_slopes_warming[1]], vert=False, positions=[2], widths=0.6, patch_artist=True,
# #                         boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
# #                         medianprops=dict(color='darkorange', linewidth=1.5),
# #                         whiskerprops=dict(color='orange'),
# #                         capprops=dict(color='orange'),
# #                         whis=(1, 99),
# #                         flierprops=dict(alpha=0.1, markeredgecolor='orange', markerfacecolor='orange', markersize=5),
# #                         showmeans=True,
# #                         meanline=True,
# #                         meanprops=dict(color='black', linestyle='--', linewidth=1.5)) 
# #         axes[1].set_yticks([])
# #         axes[1].set_xlabel('warming rate [K d$^{-1}$]')
# #         if airstream in ['BCAR', 'nonBCAR']:
# #             axes[1].set_xlim(0, 80)
# #             axes[1].set_xticks(np.arange(0, 81, 10))
# #         axes[1].grid(True, axis='x', alpha=0.3)
# #         plt.tight_layout()
# #         plt.savefig(f"{figpath}hist_box_warming_slope_d_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# # # %%
# # # For each airstream, calculate the mean and standard deviation of the slopes
# # for region in regions:
# #     print('--------------------------------------------')
# #     print('--------------------------------------------')
# #     print(f"Region: {region}")
# #     for timeframe in timeframes:
# #         print('--------------------------------------------')
# #         print(f"Timeframe: {timeframe}")
# #         path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
# #         df = pd.read_csv(path)
# #         for airstream in airstreams:
# #             df_air = df[df['airstream'] == airstream]
# #             if df_air.empty:
# #                 continue
# #             mean_cooling = df_air['cooling'].mean()
# #             std_cooling = df_air['cooling'].std()
# #             mean_warming = df_air['warming'].mean()
# #             std_warming = df_air['warming'].std()
# #             print(f"{airstream}: Cooling: {mean_cooling:.3f} ± {std_cooling:.3f} K, Warming: {mean_warming:.3f} ± {std_warming:.3f} K")
# #             mean_cooling_time = df_air['cooling_time'].mean()
# #             std_cooling_time = df_air['cooling_time'].std()
# #             mean_warming_time = df_air['warming_time'].mean()
# #             std_warming_time = df_air['warming_time'].std()
# #             print(f"{airstream}: Cooling time: {mean_cooling_time:.3f} ± {std_cooling_time:.3f} h, Warming time: {mean_warming_time:.3f} ± {std_warming_time:.3f} h")
# #             mean_cooling_slope = df_air['slope_cooling'].mean()
# #             std_cooling_slope = df_air['slope_cooling'].std()
# #             mean_warming_slope = df_air['slope_warming'].mean()
# #             std_warming_slope = df_air['slope_warming'].std()
# #             print(f"{airstream}: Cooling slope: {mean_cooling_slope:.3f} ± {std_cooling_slope:.3f} K/h, Warming slope: {mean_warming_slope:.3f} ± {std_warming_slope:.3f} K/h")


# # #%%
# # # Save summary statistics to CSV for each region and timeframe
# # for region in regions:
# #     output_rows = []
# #     for timeframe in timeframes:
# #         path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
# #         df = pd.read_csv(path)
# #         for airstream in airstreams:
# #             df_air = df[df['airstream'] == airstream]
# #             if df_air.empty:
# #                 continue
# #             mean_cooling = df_air['cooling'].mean()
# #             std_cooling = df_air['cooling'].std()
# #             mean_warming = df_air['warming'].mean()
# #             std_warming = df_air['warming'].std()
# #             mean_cooling_time = df_air['cooling_time'].mean()
# #             std_cooling_time = df_air['cooling_time'].std()
# #             mean_warming_time = df_air['warming_time'].mean()
# #             std_warming_time = df_air['warming_time'].std()
# #             mean_cooling_slope = df_air['slope_cooling'].mean()
# #             std_cooling_slope = df_air['slope_cooling'].std()
# #             mean_warming_slope = df_air['slope_warming'].mean()
# #             std_warming_slope = df_air['slope_warming'].std()
# #             output_rows.append({
# #                 'timeframe': timeframe,
# #                 'airstream': airstream,
# #                 'mean_cooling': mean_cooling,
# #                 'std_cooling': std_cooling,
# #                 'mean_warming': mean_warming,
# #                 'std_warming': std_warming,
# #                 'mean_cooling_time': mean_cooling_time,
# #                 'std_cooling_time': std_cooling_time,
# #                 'mean_warming_time': mean_warming_time,
# #                 'std_warming_time': std_warming_time,
# #                 'mean_cooling_slope': mean_cooling_slope,
# #                 'std_cooling_slope': std_cooling_slope,
# #                 'mean_warming_slope': mean_warming_slope,
# #                 'std_warming_slope': std_warming_slope
# #             })
# #     # Write to CSV
# #     out_csv = f"{datapath}TH_stats_summary_{region}.csv"
# #     pd.DataFrame(output_rows).to_csv(out_csv, index=False)

# # #%%
# # # Cleaner version of the above code, summarizing statistics for each region and timeframe
# # for region in regions:
# #     output_rows = []
# #     for airstream in airstreams:
# #         for timeframe in timeframes:
# #             path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
# #             df = pd.read_csv(path)
# #             df_air = df[df['airstream'] == airstream]
# #             if df_air.empty:
# #                 continue
# #             mean_cooling = df_air['cooling'].mean()
# #             std_cooling = df_air['cooling'].std()
# #             mean_warming = df_air['warming'].mean()
# #             std_warming = df_air['warming'].std()
# #             mean_cooling_time = df_air['cooling_time'].mean()
# #             std_cooling_time = df_air['cooling_time'].std()
# #             mean_warming_time = df_air['warming_time'].mean()
# #             std_warming_time = df_air['warming_time'].std()
# #             mean_cooling_slope = df_air['slope_cooling'].mean()
# #             std_cooling_slope = df_air['slope_cooling'].std()
# #             mean_warming_slope = df_air['slope_warming'].mean()
# #             std_warming_slope = df_air['slope_warming'].std()
# #             output_rows.append({
# #                 'timeframe': timeframe,
# #                 'airstream': airstream,
# #                 'cooling [K]': f'{mean_cooling:.2f} ± {std_cooling:.2f}',
# #                 'warming [K]': f'{mean_warming:.2f} ± {std_warming:.2f}',
# #                 'cooling time [h]': f'{mean_cooling_time:.2f} ± {std_cooling_time:.2f}',
# #                 'warming time [h]': f'{mean_warming_time:.2f} ± {std_warming_time:.2f}',
# #                 'cooling slope [K/h]': f'{mean_cooling_slope:.2f} ± {std_cooling_slope:.2f}',
# #                 'warming slope [K/h]': f'{mean_warming_slope:.2f} ± {std_warming_slope:.2f}'
# #             })
# #     # Write to CSV
# #     out_csv = f"{datapath}TH_stats_summary_{region}.csv"
# #     pd.DataFrame(output_rows).to_csv(out_csv, index=False)
# # # %%
# # # For each airsteam, plot the distribution of present and future
# # for region in regions:
# #     print('--------------------------------------------')
# #     print('--------------------------------------------')
# #     print(f"Region: {region}")
# #     for airstream in airstreams:
# #         print('--------------------------------------------')
# #         print(f"Airstream: {airstream}")
# #         dfs = {}
# #         for timeframe in timeframes:
# #             path = f"{datapath}TH_summary_{region}_{timeframe}.csv"
# #             dfs[timeframe] = pd.read_csv(path)
# #         df_present = dfs['present'][dfs['present']['airstream'] == airstream]
# #         df_future = dfs['future'][dfs['future']['airstream'] == airstream]

# #         plt.figure(figsize=(10, 6))
# #         plt.hist(df_present['cooling'], bins=50, alpha=0.5, label='present Cooling', color='blue')
# #         plt.hist(df_future['cooling'], bins=50, alpha=0.5, label='Future Cooling', color='orange')
# #         plt.xlabel('Cooling Temperature Change (K)')
# #         plt.ylabel('frequency')
# #         plt.title(f'{airstream} - {region} - Cooling Distribution')
# #         plt.legend()
# #         # plt.savefig(f"{figpath}hist_cooling_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()

# #         plt.figure(figsize=(10, 6))
# #         plt.hist(df_present['warming'], bins=50, alpha=0.5, label='present Warming', color='blue')
# #         plt.hist(df_future['warming'], bins=50, alpha=0.5, label='Future Warming', color='orange')
# #         plt.xlabel('Warming Temperature Change (K)')
# #         plt.ylabel('frequency')
# #         plt.title(f'{airstream} - {region} - Warming Distribution')
# #         plt.legend()
# #         # plt.savefig(f"{figpath}hist_warming_{airstream}_{region}.png", dpi=300, bbox_inches='tight')
# #         plt.show()
# # # %%
# # def plot_change_comparison(region):
# #     """Plot changes in cooling/warming metrics between present and future"""
# #     path = f"{datapath}TH_stats_summary_{region}.csv"
# #     df = pd.read_csv(path)
    
# #     # Extract numeric values from formatted strings
# #     def extract_mean(formatted_str):
# #         """Extract mean value from '1.23 ± 0.45' format"""
# #         return float(formatted_str.split(' ± ')[0])
    
# #     metrics = ['cooling [K]', 'warming [K]', 'cooling time [h]', 'warming time [h]', 
# #                'cooling slope [K/h]', 'warming slope [K/h]']
# #     metric_labels = ['Cooling [K]', 'Warming [K]', 'Cooling Time [h]', 'Warming Time [h]',
# #                     'Cooling Rate [K/h]', 'Warming Rate [K/h]']
    
# #     fig, axes = plt.subplots(2, 3, figsize=(18, 12))
# #     axes = axes.flatten()
    
# #     for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
# #         ax = axes[i]
        
# #         # Calculate changes for each airstream
# #         changes = []
# #         airstream_names = []
        
# #         for airstream in airstreams:
# #             present_data = df[(df['airstream'] == airstream) & (df['timeframe'] == 'present')]
# #             future_data = df[(df['airstream'] == airstream) & (df['timeframe'] == 'future')]
            
# #             if len(present_data) > 0 and len(future_data) > 0:
# #                 present_val = extract_mean(present_data[metric].values[0])
# #                 future_val = extract_mean(future_data[metric].values[0])
# #                 change = future_val - present_val
# #                 changes.append(change)
# #                 airstream_names.append(airstream)
        
# #         # Color bars based on positive/negative change
# #         colors = ['red' if x > 0 else 'blue' for x in changes]
# #         bars = ax.bar(airstream_names, changes, color=colors, alpha=0.7)
        
# #         ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
# #         ax.set_ylabel(f'Change in {label}', fontsize=12)
# #         ax.set_title(f'{region}: Future - Present', fontsize=14)
# #         ax.tick_params(axis='x', rotation=45)
        
# #         # Add value labels on bars
# #         for bar, change in zip(bars, changes):
# #             height = bar.get_height()
# #             ax.text(bar.get_x() + bar.get_width()/2., height,
# #                    f'{change:.3f}', ha='center', va='bottom' if height > 0 else 'top')
    
# #     plt.tight_layout()
# #     plt.savefig(f"{figpath}change_comparison_{region}.png", dpi=300, bbox_inches='tight')
# #     plt.show()

# # def plot_before_after_comparison(region):
# #     """Side-by-side comparison of present vs future with error bars"""
# #     path = f"{datapath}TH_stats_summary_{region}.csv"
# #     df = pd.read_csv(path)
    
# #     def extract_mean_std(formatted_str):
# #         """Extract mean and std from '1.23 ± 0.45' format"""
# #         parts = formatted_str.split(' ± ')
# #         return float(parts[0]), float(parts[1])
    
# #     fig, axes = plt.subplots(2, 3, figsize=(18, 12))
# #     axes = axes.flatten()
    
# #     metrics = ['cooling [K]', 'warming [K]', 'cooling time [h]', 
# #                'warming time [h]', 'cooling slope [K/h]', 'warming slope [K/h]']
    
# #     for idx, metric in enumerate(metrics):
# #         ax = axes[idx]
        
# #         x = np.arange(len(airstreams))
# #         width = 0.35
        
# #         present_means = []
# #         present_stds = []
# #         future_means = []
# #         future_stds = []
# #         labels = []
        
# #         for airstream in airstreams:
# #             present_data = df[(df['airstream'] == airstream) & (df['timeframe'] == 'present')]
# #             future_data = df[(df['airstream'] == airstream) & (df['timeframe'] == 'future')]
            
# #             if len(present_data) > 0 and len(future_data) > 0:
# #                 p_mean, p_std = extract_mean_std(present_data[metric].values[0])
# #                 f_mean, f_std = extract_mean_std(future_data[metric].values[0])
                
# #                 present_means.append(p_mean)
# #                 present_stds.append(p_std)
# #                 future_means.append(f_mean)
# #                 future_stds.append(f_std)
# #                 labels.append(airstream)
        
# #         x_adj = np.arange(len(labels))
# #         ax.bar(x_adj - width/2, present_means, width, yerr=present_stds, 
# #                label='Present', color='blue', alpha=0.7, capsize=5)
# #         ax.bar(x_adj + width/2, future_means, width, yerr=future_stds,
# #                label='Future', color='orange', alpha=0.7, capsize=5)
        
# #         ax.set_ylabel(metric, fontsize=12)
# #         ax.set_title(f'{region}: {metric}')
# #         ax.set_xticks(x_adj)
# #         ax.set_xticklabels(labels)
# #         ax.legend()
# #         ax.grid(True, alpha=0.3)
    
# #     plt.tight_layout()
# #     plt.savefig(f"{figpath}before_after_{region}.png", dpi=300, bbox_inches='tight')
# #     plt.show()

# # def plot_slope_visualization(region):
# #     """Show cooling/warming as slopes on a temperature-time diagram"""
# #     path = f"{datapath}TH_stats_summary_{region}.csv"
# #     df = pd.read_csv(path)
    
# #     def extract_mean(formatted_str):
# #         return float(formatted_str.split(' ± ')[0])
    
# #     fig, axes = plt.subplots(2, 2, figsize=(16, 10))
# #     axes = axes.flatten()
    
# #     available_airstreams = df['airstream'].unique()
    
# #     for i, airstream in enumerate(available_airstreams[:4]):  # Limit to 4 for subplot layout
# #         ax = axes[i]
        
# #         for timeframe, color in colors.items():
# #             data = df[(df['airstream'] == airstream) & (df['timeframe'] == timeframe)]
# #             if len(data) == 0: 
# #                 continue
            
# #             # Extract values using the correct column names
# #             cooling_time = extract_mean(data['cooling time [h]'].values[0])
# #             warming_time = extract_mean(data['warming time [h]'].values[0])
# #             cooling_slope = extract_mean(data['cooling slope [K/h]'].values[0])
# #             warming_slope = extract_mean(data['warming slope [K/h]'].values[0])
# #             cooling_amount = extract_mean(data['cooling [K]'].values[0])
# #             warming_amount = extract_mean(data['warming [K]'].values[0])
            
# #             # Time points (cooling_time is from start to minimum)
# #             t1 = -168  # Start
# #             t2 = -cooling_time  # Time of minimum TH
# #             t3 = 0  # End
            
# #             # Temperature points (using actual cooling/warming amounts)
# #             TH1 = 0  # Start temperature (relative)
# #             TH2 = cooling_amount  # Change during cooling (negative)
# #             TH3 = TH2 + warming_amount  # Final temperature
            
# #             # Plot the trajectory
# #             times = [t1, t2, t3]
# #             temps = [TH1, TH2, TH3]
            
# #             ax.plot(times, temps, 'o-', color=color, linewidth=3, 
# #                    label=f'{timeframe}', markersize=8)
            
# #             # Add slope annotations
# #             ax.annotate(f'{cooling_slope:.3f} K/h', 
# #                        xy=((t1+t2)/2, (TH1+TH2)/2), 
# #                        xytext=(10, 10), textcoords='offset points',
# #                        fontsize=10, color=color)
# #             ax.annotate(f'{warming_slope:.3f} K/h', 
# #                        xy=((t2+t3)/2, (TH2+TH3)/2),
# #                        xytext=(10, 10), textcoords='offset points', 
# #                        fontsize=10, color=color)
        
# #         ax.set_xlabel('Time [h]', fontsize=12)
# #         ax.set_ylabel('Temperature Change [K]', fontsize=12)
# #         ax.set_title(f'{airstream} - {region}', fontsize=14)
# #         ax.legend()
# #         ax.grid(True, alpha=0.3)
# #         ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    
# #     # Hide unused subplots
# #     for j in range(len(available_airstreams), 4):
# #         axes[j].set_visible(False)
    
# #     plt.tight_layout()
# #     plt.savefig(f"{figpath}slope_visualization_{region}.png", dpi=300, bbox_inches='tight')
# #     plt.show()

# # # Usage:
# # for region in regions:
# #     plot_change_comparison(region)
# #     plot_before_after_comparison(region)
# #     plot_slope_visualization(region)


# # # %%
# # # extra stuff
# # def calculate_trajectory_based_stats(df, varname):
# #     """Calculate stats using the same trajectory-based approach as your minima analysis"""
# #     # First, get the same trajectories you analyze individually
# #     valid_trajs = []
# #     traj_ids = []
    
# #     for traj_id, group in df.groupby('traj_id'):
# #         # Apply the same filtering criteria you use in your minima analysis
# #         if len(group) > 0 and not group[varname].isna().all():
# #             valid_trajs.append(group)
# #             traj_ids.append(traj_id)
    
# #     # Now calculate ensemble stats using only these valid trajectories
# #     if not valid_trajs:
# #         return [], [], [], []
    
# #     # Combine all valid trajectory data
# #     combined_df = pd.concat(valid_trajs, ignore_index=True)
    
# #     # Calculate stats using the same approach as before
# #     stats = {}
# #     for t, group in combined_df.groupby('rel_time'):
# #         vals = group[varname]
# #         stats[t] = (
# #             vals.mean(),
# #             vals.quantile(0.25), 
# #             vals.quantile(0.75)
# #         )
    
# #     # Sort by rel_time
# #     rel_times = sorted(stats.keys())
# #     means = [stats[t][0] for t in rel_times]
# #     q25s = [stats[t][1] for t in rel_times]
# #     q75s = [stats[t][2] for t in rel_times]
    
# #     return rel_times, means, q25s, q75s, traj_ids

# # def analyze_trajectory_minima_with_ensemble(df, varname='TH'):
# #     """Combined analysis that ensures consistency between ensemble and individual trajectory stats"""
    
# #     # Get ensemble stats and valid trajectory IDs
# #     rel_times, means, q25s, q75s, valid_traj_ids = calculate_trajectory_based_stats(df, varname)
    
# #     # Analyze individual trajectories (using same trajectories as ensemble)
# #     summary_list = []
# #     for traj_id in valid_traj_ids:
# #         group = df[df['traj_id'] == traj_id]
        
# #         th_0 = group.loc[np.isclose(group['rel_time'], 0.0), varname]
# #         th_168 = group.loc[np.isclose(group['rel_time'], -168.0), varname]
# #         th_min = group[varname].min()
# #         th_min_time = group.loc[group[varname] == th_min, 'rel_time'].values[0] if not group[group[varname] == th_min].empty else np.nan
        
# #         summary_list.append({
# #             'traj_id': traj_id,
# #             f'{varname}_at_0': th_0.iloc[0] if not th_0.empty else np.nan,
# #             f'{varname}_at_-168': th_168.iloc[0] if not th_168.empty else np.nan,
# #             f'{varname}_min': th_min,
# #             f'{varname}_min_time': th_min_time
# #         })
    
# #     individual_stats = pd.DataFrame(summary_list)
    
# #     return {
# #         'ensemble': (rel_times, means, q25s, q75s),
# #         'individual': individual_stats,
# #         'n_trajectories': len(valid_traj_ids)
# #     }
# %%
# # --- First pass: determine auto y-limits for each airstream using trajectory-based stats ---
# auto_limits = {airstream: [] for airstream in airstreams}
# varname = 'TH'
# rounding = 5

# for region in regions:
#     dfs = {}
#     for timeframe in timeframes:
#         path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
#         dfs[timeframe] = pd.read_csv(path)
    
#     for airstream in airstreams:
#         fig, ax = plt.subplots()
#         for timeframe in timeframes:
#             df = dfs[timeframe]
#             df_air = df[df['airstream'] == airstream]
#             if df_air.empty:
#                 continue
            
#             # Use the combined analysis to get consistent trajectory sets
#             results = analyze_trajectory_minima_with_ensemble(df_air, 'TH')
#             rel_times, means, q25s, q75s = results['ensemble']
            
#             ax.plot(rel_times, means, color=colors[timeframe])
#             ax.fill_between(rel_times, q25s, q75s, alpha=0.2, color=colors[timeframe])
        
#         # Record the auto y-limits
#         auto_limits[airstream].append(ax.get_ylim())
#         plt.close(fig)

# # Find the global min/max of the auto-limits for each airstream, round to 5
# th_limits = {}
# for airstream in airstreams:
#     if auto_limits[airstream]:  # Check if we have data for this airstream
#         ymins, ymaxs = zip(*auto_limits[airstream])
#         th_min = rounding * np.floor(min(ymins) / rounding)
#         th_max = rounding * np.ceil(max(ymaxs) / rounding)
#         th_limits[airstream] = (th_min, th_max)
#     else:
#         # Default limits if no data available
#         th_limits[airstream] = (250, 350)

# # --- Second pass: plot with consistent y-limits using trajectory-based analysis ---
# for region in regions:
#     dfs = {}
#     for timeframe in timeframes:
#         path = f"/net/litho/atmosdyn2/nfieldhouse/data/lagranto/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
#         dfs[timeframe] = pd.read_csv(path)
    
#     for airstream in airstreams:
#         stats_dict = {}
#         minima_dict = {}
        
#         for timeframe in timeframes:
#             df = dfs[timeframe]
#             df_air = df[df['airstream'] == airstream]
#             if df_air.empty:
#                 continue
            
#             # Use the combined analysis
#             results = analyze_trajectory_minima_with_ensemble(df_air, 'TH')
#             stats_dict[timeframe] = results['ensemble']
#             minima_dict[timeframe] = results['individual']
            
#             print(f"{region} {timeframe} {airstream}: {results['n_trajectories']} trajectories analyzed")
        
#         if stats_dict:
#             # Plot ensemble evolution with consistent y-limits
#             plot_weighted_stats(stats_dict, th_limits, region, airstream, 'TH')
            
#             # Now your minima analysis uses the exact same trajectories
#             for timeframe in timeframes:
#                 if timeframe in minima_dict:
#                     individual_stats = minima_dict[timeframe]
#                     avg_min_time = individual_stats['TH_min_time'].mean()
#                     avg_min_TH = individual_stats['TH_min'].mean()
                    
#                     print(f"{region} {timeframe} {airstream}:")
#                     print(f"  Average minimum time: {avg_min_time:.1f} hours")
#                     print(f"  Average minimum TH: {avg_min_TH:.2f} K")
# %%

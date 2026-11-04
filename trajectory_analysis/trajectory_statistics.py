"""Evolution of the potential temperature along back-trajectories of marine cold-air outbreaks.

Plots the weighted mean evolution of the potential temperature (TH) along the trajectories and
the distribution of the cooling and warming rates, for the present and the future climate.
The input files are written by trajectory_statistics_writesummary.py and
trajectories/write_trajectories.py (see the README).

Usage (from any directory):

    python trajectory_statistics.py FIGURE [FIGURE ...]

where FIGURE is one of: F3, S4, S5 (see FIGURES below). S4 also prints the numbers of the table
in Fig. S4c, f (cooling and warming counts, shares and mean rates). The input data are read from
$CAO_DATA_ROOT/lagranto and the figures are written to $CAO_PLOTS_ROOT/traj
(see config/config.py).
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
from scipy.stats import ttest_ind
import sys
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
import config

plt.rcParams.update({
    'font.size': 16,           # Default font size for ticks
    'axes.titlesize': 20,      # Title font size
    'axes.labelsize': 20,      # X and Y label font size
    'legend.fontsize': 20,     # Legend font size
    'xtick.labelsize': 16,     # X-axis tick font size
    'ytick.labelsize': 16      # Y-axis tick font size
})
plt.rcParams["font.sans-serif"] = "Arial"

pt = 1./72.27 # Hundreds of years of history... 72.27 points to an inch.

textwidth = 418.25368

figsizes = {"onecol": textwidth*pt, "twocol": 0.48*textwidth*pt, "threecol":0.3*textwidth*pt, "8":8, "10":10}

my_width = figsizes["twocol"]

fourto3 = 4/3

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

def plot_combined_evolution_rates_cooling_warming(region, airstream, varname='TH'):
    """
    Create a 1x3 plot combining:
    1. Variable evolution (like plot_weighted_stats)
    2. Vertical boxplot for cooling rates
    3. Vertical boxplot for warming rates
    """

    arial_user_font = os.path.expanduser('~/.local/share/fonts/truetype/arial/arial.ttf')
    if os.path.exists(arial_user_font):
        try:
            matplotlib.font_manager.fontManager.addfont(arial_user_font)
        except Exception:
            pass
    available_fonts = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
    preferred_font = 'Arial' if 'Arial' in available_fonts else 'Liberation Sans'

    plt.rcParams.update({
    'font.size': 6,           # Default font size for ticks
    'axes.titlesize': 6,      # Title font size
    'axes.labelsize': 6,      # X and Y label font size
    'legend.fontsize': 6,     # Legend font size
    'xtick.labelsize': 4,     # X-axis tick font size
    'ytick.labelsize': 4,     # Y-axis tick font size
    'font.family': 'sans-serif',
    'font.sans-serif': [preferred_font, 'Arial', 'Liberation Sans', 'DejaVu Sans']
})
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
        ax1.plot(rel_times, means, label=f'mean {varname} ({timeframe})' if varname != 'TH' else f'mean Θ ({timeframe})', color=colors[timeframe], linewidth=0.8)
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
        ax1.axvline(x=min_th_time, color=colors[timeframe], linestyle='--' if timeframe == 'future' else '-', alpha=0.7, linewidth=0.8)
    
    ax1.set_xlabel(r'days relative to t$_{\mathrm{MCAO}}$')
    if varname == 'Q':
        ax1.set_ylabel(f'q [g kg$^{{-1}}$]')
        ax1.set_ylim(0, 5)
    elif varname == 'TH':
        ax1.set_ylabel('Θ [K]')
        ax1.set_ylim(250, 300)
    else:
        ax1.set_ylabel(f'{varname} [{units_dict[varname]}]')
    
    ax1.set_xlim(-168, 0)
    ax1.set_xticks(np.arange(-168, 1, 24))
    ax1.set_xticklabels([str(int(x//24)) for x in np.arange(-168, 1, 24)])
    ax1.set_xticks(np.arange(-168, 1, 6), minor=True)
    ax1.grid(True, alpha=0.3)
    legend = ax1.legend()
    # ax1.set_title(f'a) {region} - {airstream}', loc='left')
    
    # Second subplot: Cooling rates boxplot
    ax2 = axes[1]
    if varname == 'TH':
        cooling_rate_present = df_air_present['cooling_rate_6h'].dropna()
        cooling_rate_future = df_air_future['cooling_rate_6h'].dropna()
        rate_label = 'Θ rate [K 6h$^{-1}$]'
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

    # Force font on text artists in this figure.
    for ax in (ax1, ax2, ax5):
        ax.xaxis.label.set_fontname(preferred_font)
        ax.yaxis.label.set_fontname(preferred_font)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname(preferred_font)
    if legend is not None:
        for text in legend.get_texts():
            text.set_fontname(preferred_font)

    # # Debug: print actual resolved font files used by representative text artists.
    # y1_props = ax1.yaxis.label.get_fontproperties()
    # y2_props = ax2.yaxis.label.get_fontproperties()
    # print(f"[FONT DEBUG] preferred_font={preferred_font}")
    # print(f"[FONT DEBUG] ax1 ylabel family={y1_props.get_family()} name={y1_props.get_name()}")
    # print(f"[FONT DEBUG] ax1 ylabel file={matplotlib.font_manager.findfont(y1_props, fallback_to_default=True)}")
    # print(f"[FONT DEBUG] ax2 ylabel family={y2_props.get_family()} name={y2_props.get_name()}")
    # print(f"[FONT DEBUG] ax2 ylabel file={matplotlib.font_manager.findfont(y2_props, fallback_to_default=True)}")
    # ticklabels = ax1.get_xticklabels()
    # if ticklabels:
    #     tick_props = ticklabels[0].get_fontproperties()
    #     print(f"[FONT DEBUG] ax1 xtick[0] family={tick_props.get_family()} name={tick_props.get_name()}")
    #     print(f"[FONT DEBUG] ax1 xtick[0] file={matplotlib.font_manager.findfont(tick_props, fallback_to_default=True)}")
    # if legend is not None and legend.get_texts():
    #     leg_props = legend.get_texts()[0].get_fontproperties()
    #     print(f"[FONT DEBUG] legend[0] family={leg_props.get_family()} name={leg_props.get_name()}")
    #     print(f"[FONT DEBUG] legend[0] file={matplotlib.font_manager.findfont(leg_props, fallback_to_default=True)}")

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
    my_width = 6.69291

    if region == 'JS':
        letters = ['a', 'b', 'c']
    elif region == 'NA':
        letters = ['a', 'b', 'c']
    
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
                    medianprops=dict(color='blue', linewidth=0),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box2 = ax1.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0),
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
                    medianprops=dict(color='blue', linewidth=0),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))

    box2 = ax2.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0),
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
                    medianprops=dict(color='blue', linewidth=0),
                    whiskerprops=dict(color='blue'),
                    capprops=dict(color='blue'),
                    whis=(1, 99),
                    flierprops=dict(alpha=0.1, markeredgecolor='blue', markerfacecolor='blue', markersize=1),
                    showmeans=True,
                    meanline=True,
                    meanprops=dict(color='black', linestyle='--', linewidth=0.8))
    box6 = ax5.boxplot([data_warming[1]], positions=[2], widths=0.6, patch_artist=True,
                    boxprops=dict(facecolor='orange', color='orange', alpha=0.5),
                    medianprops=dict(color='darkorange', linewidth=0),
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
    path = f"{config.LAGRANTO_DIR}/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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
        path = f"{config.LAGRANTO_DIR}/CAO_ALL_TRAJS_{region}_{timeframe}_ALL.csv"
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


def figure_s4():
    """Distribution of the cooling rates and evolution of TH, NA low-level air stream, present and
    future (Fig. S4a, b, d, e). The table (Fig. S4c, f) is printed by cooling_rate_distribution.
    The figures use the default matplotlib font sizes."""
    plt.rcParams.update(matplotlib.rcParamsDefault)
    plt.rcParams["font.sans-serif"] = "Arial"
    for timeframe in timeframes:
        compare_cooling_rates_fixed_vs_ensemble("NA", "LL", timeframe, -96, -90)
    cooling_rate_distribution("NA", "LL", -96, -90)


# Figure name -> function that produces it
FIGURES = {
    # Evolution of TH and cooling/warming rates of the low-level air stream, Japan Sea and northwest Atlantic
    "F3": lambda: [plot_combined_evolution_rates_cooling_warming(region, "LL") for region in ("JS", "NA")],
    "S4": figure_s4,
    # Warming rates split into WBC, non-WBC and land (three panels), northwest Atlantic
    "S5": lambda: plot_combined_evolution_rates_wbc_nonwbc_land("NA", "LL"),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("figures", nargs="+", choices=sorted(FIGURES), help="figures to produce")
    args = parser.parse_args()
    os.makedirs(figpath, exist_ok=True)
    for name in args.figures:
        print(f"Producing {name}")
        with plt.rc_context():  # the plot settings of one figure must not affect the next
            FIGURES[name]()


if __name__ == "__main__":
    main()

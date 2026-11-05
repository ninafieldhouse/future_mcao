# Marine cold-air outbreak (MCAO) analysis scripts

Scripts for identifying marine cold-air outbreaks (MCAOs) in a CESM2 large ensemble, tracing the air masses that feed them with Lagrangian back-trajectories, and comparing the present (1990-2000) and future (2090-2100, SSP3-7.0) climates. Two regions are analysed: the Japan-adjacent seas (`JS`, 130-150E, 30-60N) and the northwest Atlantic (`NA`, 75-65W, 30-45N). The ensemble has five members (`1000`, `1100`, `1200`, `1300`, `1400`). Only winter (DJF) is analysed.

## Layout

| Directory | Purpose |
|---|---|
| `config/` | Machine-specific paths (`config.sh`, `config.py`) |
| `climatology/` | Cold air outbreak index (CAOI), event detection, event statistics and figures |
| `trajectories/` | Trajectory start points, Lagranto runs, trajectory post-processing |
| `trajectory_analysis/` | Trajectory density, plots and statistics |

## Requirements

- [CESM2](https://www.cesm.ucar.edu/models/cesm2/download) output on the CESM archive (CAM `h3` and CICE `h1` history files)
- [CDO](https://code.mpimet.mpg.de/projects/cdo)
- [LAGRANTO](https://iacweb.ethz.ch/staff/sprenger/lagranto/) with the CESM interface (`caltra`, `density`)
- Python 3 with numpy, pandas, xarray, netCDF4, scipy, matplotlib, cartopy, cmcrameri, metpy, geopy, numba, `dypy` and the `lagranto` Python package. The figure scripts need a recent environment: they were run with Python 3.13, matplotlib 3.10, numpy 2.3, xarray 2025.10, cartopy 0.25 and cmcrameri 1.8 (matplotlib 3.5 or newer is required for `layout='constrained'`). The figure scripts do not need `lagranto` or `dypy`.

## Configuration

Machine-specific paths are set in `config/config.sh` (for shell scripts) and `config/config.py` (for Python scripts). Each can be overridden by exporting the variable before running a script.

| Variable | Meaning |
|---|---|
| `CAO_DATA_ROOT` | Root of all derived data (per-member subdirectories `<member_id>/...`) |
| `CAO_CESM_ARCHIVE` | CESM2 archive containing the raw model output |
| `CAO_RUN_DIR` | Scratch directory for CDO and Lagranto runs |
| `CAO_PLOTS_ROOT` | Output directory for figures |

Scientific parameters (members, regions, timeframes, thresholds, region boxes) are written as literals in each script and are not part of the configuration.

Most scripts take the ensemble member as their first argument, for example `bash climatology/CAO_index.sh 1000`.

## Workflow

### Branch A: MCAO climatology and statistics (`climatology/`)

Run once per ensemble member, in this order.

1. `CAO_index.sh` computes the CAOI from CESM output for each region and timeframe, using CDO.
2. `CAO_index_clim.sh` computes the six-hourly climatology and the anomalies of the CAOI.
3. `CAO_index_evolution.py` averages the anomalies over the region box (latitude weighted), writes the CAOI time series `CAOI_<timeframe>.csv` and plots it.
4. `extract_events_CAOI_6hourly.sh` detects the events. Its final output, `summary_events_DJF90_<timeframe>_overstep_24_CAOI_<region>_filtered`, is the event list used by the trajectory branch.
5. `CAO_index_statistics_clean_write_data.py` and `CAO_index_statistics_clean_write_precip_data.py` write the event-based fields (CAOI and precipitation) to files.
6. `CAO_index_statistics_clean_mulitimember.py` combines the five members and produces the ensemble figures (see [Paper figures](#paper-figures)).

### Branch B: trajectories (`trajectories/`, `trajectory_analysis/`)

Run `trajectories/` steps per ensemble member, except where noted.

1. `CAO_masks_edit.sh` calls `CAO_masks_2.py` for every event. It writes the MCAO mask and the Lagranto start file for the event date.
2. `calculate_trajs.sh` runs the four `loop_through_events_<region>_<timeframe>.sh` scripts in parallel. Each loops over all events and calls the matching `backward_forward_trajectories_cao_<region>_<timeframe>.sh`, which runs Lagranto `caltra`.
   `backward_forward_trajectories_cao_singleevent.sh` runs one event only.
3. `find_missing_trajs.sh` lists events that have no trajectory file, so that they can be rerun.
4. `add_th_column.sh` adds potential temperature `TH` to every trajectory file (`B_*_with_TH.txt`).
5. `write_trajectories.py` reads all five members, classifies the trajectories into air streams, and writes one concatenated file per air stream, region and timeframe (`CAO_<airstream>_trajs_<region>_<timeframe>_ALL.txt`).
6. `fix_trajectory_times.py` makes the time axis consistent (29 steps, 0 to -168 h in steps of -6 h).
7. `add_dth_column.py` adds the potential-temperature tendencies `dTH` and `dTH_24`.
8. `combine_traj.sh` is a helper that concatenates trajectory files and keeps only the first header.

Then, in `trajectory_analysis/`:

- `density_trajs.sh` computes trajectory densities on a grid with Lagranto `density`. The Lagranto `density` tool only accepts file names without a directory part, so the script runs inside `$CAO_LAGRANTO_DIR` and then moves the results to `$CAO_DENSITY_DIR`.
- `plot_density.py` plots the 24-hourly potential temperature change of the low-level air stream and its present-future difference (see [Paper figures](#paper-figures)).
- `plot_trajectories.py` plots trajectories by air stream and the air-stream fractions.
- `trajectory_statistics_writesummary.py` writes a summary file per region and timeframe (`TH_summary_<region>_<timeframe>.csv`) that `trajectory_statistics.py` reads.
- `trajectory_statistics.py` plots the evolution of the potential temperature along the trajectories and the distribution of the cooling and warming rates, with significance tests (see [Paper figures](#paper-figures)).

## Paper figures

The figure scripts read the data from `$CAO_DATA_ROOT` and write the figures to `$CAO_PLOTS_ROOT`.

| Figure | Command |
|---|---|
| Fig. 1 | `python climatology/CAO_index_statistics_clean_mulitimember.py F1` |
| Fig. S1 | `python climatology/CAO_index_statistics_clean_mulitimember.py S1` |
| Fig. S2 | `python climatology/CAO_index_statistics_clean_mulitimember.py S2` |
| Fig. 2 | `python climatology/CAO_index_statistics_clean_mulitimember.py F2` |
| Fig. S3 | `python climatology/CAO_index_statistics_clean_mulitimember.py S3` |
| Fig. 3 | `python trajectory_analysis/trajectory_statistics.py F3` |
| Fig. S4 | `python trajectory_analysis/trajectory_statistics.py S4` |
| Fig. S5 | `python trajectory_analysis/trajectory_statistics.py S5` |
| Fig. 4 | `python trajectory_analysis/plot_density.py F4_maps F4_diff` |

## Notes

- Region and timeframe loops are hard-coded inside the shell scripts. To process a subset, edit the `regions` and `timeframes` arrays at the top of the script.
- Archived and superseded scripts are not part of this repository.

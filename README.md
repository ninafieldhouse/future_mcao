# Marine cold-air outbreak (MCAO) analysis scripts

Scripts for identifying marine cold-air outbreaks (CAOs) in a CESM2 large ensemble, tracing the air masses that feed them with Lagrangian back-trajectories, and comparing the present (1990-2000) and future (2090-2100, SSP3-7.0) climates. Two regions are analysed: the Japan Sea (`JS`, 130-150E, 30-60N) and the northwest Atlantic (`NA`, 75-65W, 30-45N). The ensemble has five members (`1000`, `1100`, `1200`, `1300`, `1400`). Only winter (DJF) is analysed.

## Layout

| Directory | Purpose |
|---|---|
| `config/` | Machine-specific paths (`config.sh`, `config.py`) |
| `climatology/` | CAO index, event detection, event statistics and figures |
| `trajectories/` | Trajectory start points, Lagranto runs, trajectory post-processing |
| `trajectory_analysis/` | Trajectory density, plots and statistics |

## Requirements

- CESM2 output on the CESM archive (CAM `h3` and CICE `h1` history files)
- [CDO](https://code.mpimet.mpg.de/projects/cdo)
- LAGRANTO with the CESM interface (`caltra`, `density`), loaded through the `dyn_tools` module
- Python 3 with numpy, pandas, xarray, netCDF4, scipy, matplotlib, cartopy, cmcrameri, metpy, geopy, numba, `dypy` and the `lagranto` Python package

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

## Method summary

**CAO index.** The index is the difference between the potential temperature of the sea surface (SST referenced to surface pressure) and of the air at about 850 hPa: `CAOI = SSTH - TH_850`. Grid points with sea ice or land are masked out. A grid point is counted as a CAO point where `CAOI > 8 K`.

**Events.** Six-hourly CAO index values in DJF are ranked and those above the 90th percentile are kept. Consecutive time steps, and time steps separated by gaps of up to 24 h, are merged into one event, starting with the most intense. Only events lasting longer than 24 h are kept.

**Trajectories.** For each event, 7-day (168 h) backward trajectories are started from all CAO grid points in the region, from the surface up to an assumed boundary-layer top of 850 hPa. Potential temperature (`TH`), its tendency (`dTH`, `dTH_24`) and other fields are traced along each trajectory.

**Air streams.** Trajectories are classified as follows.

| Name | Definition |
|---|---|
| `LL` | Low-level: pressure above 700 hPa at all 29 time steps |
| `UL` | All other trajectories |
| `BCAR` | `LL` trajectories with `TH < 280 K` over land (or over sea ice in a box over the northwest Atlantic) for more than 12 consecutive time steps, and sea-level pressure above 1020 hPa at least once |
| `nonBCAR` | `LL` trajectories that are not `BCAR` |
| `DI` | Trajectories with a pressure change of more than 400 hPa relative to their start within the first 12 time steps |

## Workflow

### Branch A: CAO climatology and statistics (`climatology/`)

Run once per ensemble member, in this order.

1. `CAO_index.sh` computes the CAO index from CESM output for each region and timeframe, using CDO.
2. `CAO_index_clim.sh` computes the six-hourly climatology and the anomalies of the index.
3. `extract_events_CAOI_6hourly.sh` detects the events. Its final output, `summary_events_DJF90_<timeframe>_overstep_24_CAOI_<region>_filtered`, is the event list used by the trajectory branch.
4. `CAO_index_statistics_clean_write_data.py` and `CAO_index_statistics_clean_write_precip_data.py` write the event-based fields (CAO index and precipitation) to files.
5. `CAO_index_statistics_clean_mulitimember.py` combines the five members and produces the event statistics and the climatology figures.
6. `CAO_index_evolution.py` plots the evolution of the index during events.

### Branch B: trajectories (`trajectories/`, `trajectory_analysis/`)

Run `trajectories/` steps per ensemble member, except where noted.

1. `CAO_masks_edit.sh` calls `CAO_masks_2.py` for every event. It writes the CAO mask and the Lagranto start file for the event date.
2. `calculate_trajs.sh` runs the four `loop_through_events_<region>_<timeframe>.sh` scripts in parallel. Each loops over all events and calls the matching `backward_forward_trajectories_cao_<region>_<timeframe>.sh`, which runs Lagranto `caltra`.
   `backward_forward_trajectories_cao_singleevent.sh` runs one event only.
3. `find_missing_trajs.sh` lists events that have no trajectory file, so that they can be rerun.
4. `add_th_column.sh` adds potential temperature `TH` to every trajectory file (`B_*_with_TH.txt`).
5. `write_trajectories.py` reads all five members, classifies the trajectories into air streams, and writes one concatenated file per air stream, region and timeframe (`CAO_<airstream>_trajs_<region>_<timeframe>_ALL.txt`).
6. `fix_trajectory_times.py` makes the time axis consistent (29 steps, 0 to -168 h in steps of -6 h).
7. `add_dth_column.py` adds the potential-temperature tendencies `dTH` and `dTH_24`.
8. `combine_traj.sh` is a helper that concatenates trajectory files and keeps only the first header.

Then, in `trajectory_analysis/`:

- `density_trajs.sh` computes trajectory densities on a grid with Lagranto `density`. Output goes to `$CAO_DENSITY_DIR`.
- `plot_density.py` plots the densities and the present-future differences.
- `plot_trajectories.py` plots trajectories by air stream and the air-stream fractions.
- `trajectory_statistics.py` and `trajectory_statistics_writesummary.py` compute the weighted statistics and significance tests (t-test, Mann-Whitney U) and plot the warming along trajectories.

## Notes

- Region and timeframe loops are hard-coded inside the shell scripts. To process a subset, edit the `regions` and `timeframes` arrays at the top of the script.
- Archived and superseded scripts are not part of this repository.

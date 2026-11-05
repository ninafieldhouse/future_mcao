#!/bin/bash
# Shared path configuration for the shell scripts. Every script sources this file near the
# top (`source "$(dirname "$0")/../config/config.sh"`), so the paths are defined in one
# place. There are no default paths: export the four variables below before running a
# script (see the Configuration section of the README), for example in ~/.bashrc:
#
#   export CAO_DATA_ROOT=/path/to/derived/data
#   export CAO_PLOTS_ROOT=/path/to/figures
#   export CAO_CESM_ARCHIVE=/path/to/cesm2/output
#   export CAO_RUN_DIR=/path/to/scratch
#
# Scientific parameters (ensemble members, regions, timeframes, airstreams, CESM case
# names, region boxes, year ranges) are fixed by the study design and are kept as
# explicit literals in each script rather than centralized here.

: "${CAO_DATA_ROOT:?CAO_DATA_ROOT is not set (directory for all derived data)}"
: "${CAO_PLOTS_ROOT:?CAO_PLOTS_ROOT is not set (output directory for figures)}"
: "${CAO_CESM_ARCHIVE:?CAO_CESM_ARCHIVE is not set (CESM2 output archive)}"
: "${CAO_RUN_DIR:?CAO_RUN_DIR is not set (scratch directory for CDO and LAGRANTO runs)}"
export CAO_DATA_ROOT CAO_PLOTS_ROOT CAO_CESM_ARCHIVE CAO_RUN_DIR
export CAO_LAGRANTO_DIR="${CAO_DATA_ROOT}/lagranto"
export CAO_DENSITY_DIR="${CAO_LAGRANTO_DIR}/density"

# Scripts that call LAGRANTO (caltra, density, newtime) run this function. It needs
# CAO_LAGRANTO_HOME, the installation directory of LAGRANTO.CESM.
cao_setup_lagranto() {
    : "${CAO_LAGRANTO_HOME:?CAO_LAGRANTO_HOME is not set (installation directory of LAGRANTO.CESM)}"
    export LAGRANTO="$CAO_LAGRANTO_HOME"
    [ -d "$CAO_LAGRANTO_HOME/bin" ] && export PATH="$CAO_LAGRANTO_HOME/bin:$PATH"
    for tool in caltra density newtime; do
        command -v "$tool" >/dev/null 2>&1 || { echo "Error: LAGRANTO tool '$tool' not found; check CAO_LAGRANTO_HOME and PATH" >&2; return 1; }
    done
}

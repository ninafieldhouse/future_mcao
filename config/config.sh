#!/bin/bash
# Shared path configuration for the CAO pipeline shell scripts. Source this near the top
# of every script (`source "$(dirname "$0")/../config/config.sh"`) instead of hardcoding
# paths locally. Override any CAO_* value by exporting it before running.
#
# Scientific parameters (ensemble members, regions, timeframes, airstreams, CESM case
# names, region boxes, year ranges) are fixed by the study design and are kept as
# explicit literals in each script rather than centralized here.

export CAO_DATA_ROOT="${CAO_DATA_ROOT:-/net/litho/atmosdyn2/nfieldhouse/data}"
export CAO_PLOTS_ROOT="${CAO_PLOTS_ROOT:-/home/nfieldhouse/ma/plots}"
export CAO_CESM_ARCHIVE="${CAO_CESM_ARCHIVE:-/net/meso/climphys/cesm212}"
export CAO_RUN_DIR="${CAO_RUN_DIR:-/home/nfieldhouse/ma/run}"
export CAO_LAGRANTO_DIR="${CAO_DATA_ROOT}/lagranto"
export CAO_DENSITY_DIR="${CAO_LAGRANTO_DIR}/density"

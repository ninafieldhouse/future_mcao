"""Shared path configuration for the CAO pipeline. These are the only things that differ
between environments/machines; override any of them by setting the corresponding
environment variable before running a script. Scientific parameters (ensemble members,
regions, timeframes, airstreams, CESM case names, region boxes) are fixed by the study
design and are kept as explicit literals in each script rather than centralized here."""

import os

DATA_ROOT = os.environ.get("CAO_DATA_ROOT", "/net/litho/atmosdyn2/nfieldhouse/data")
PLOTS_ROOT = os.environ.get("CAO_PLOTS_ROOT", "/home/nfieldhouse/ma/plots")
CESM_ARCHIVE = os.environ.get("CAO_CESM_ARCHIVE", "/net/meso/climphys/cesm212")
RUN_DIR = os.environ.get("CAO_RUN_DIR", "/home/nfieldhouse/ma/run")

LAGRANTO_DIR = os.path.join(DATA_ROOT, "lagranto")
DENSITY_DIR = os.path.join(LAGRANTO_DIR, "density")

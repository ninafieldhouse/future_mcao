"""Shared path configuration for the Python scripts. There are no default paths: export
the four variables below before running a script (see the Configuration section of the
README). Scientific parameters (ensemble members, regions, timeframes, airstreams, CESM
case names, region boxes) are fixed by the study design and are kept as explicit literals
in each script rather than centralized here."""

import os


def _require(name, description):
    try:
        return os.environ[name]
    except KeyError:
        raise RuntimeError(
            f"{name} is not set ({description}). "
            "See the Configuration section of the README."
        ) from None


DATA_ROOT = _require("CAO_DATA_ROOT", "directory for all derived data")
PLOTS_ROOT = _require("CAO_PLOTS_ROOT", "output directory for figures")
CESM_ARCHIVE = _require("CAO_CESM_ARCHIVE", "CESM2 output archive")
RUN_DIR = _require("CAO_RUN_DIR", "scratch directory for CDO and LAGRANTO runs")

LAGRANTO_DIR = os.path.join(DATA_ROOT, "lagranto")
DENSITY_DIR = os.path.join(LAGRANTO_DIR, "density")

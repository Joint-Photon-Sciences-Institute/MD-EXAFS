"""MD-EXAFS: A Python package for calculating EXAFS spectra from MD trajectories."""

__version__ = "0.1.0"

from .config import load_config, ConfigError
from .trajectory_processor import process_trajectory
from .chi_averager import average_chi_files

__all__ = ["load_config", "ConfigError", "process_trajectory", "average_chi_files"]
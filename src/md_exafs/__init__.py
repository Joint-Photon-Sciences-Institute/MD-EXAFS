"""MD-EXAFS: A Python package for calculating EXAFS spectra from MD trajectories."""

__version__ = "2.0.0"

from .config import load_config, ConfigError
from .trajectory_processor import process_trajectory
from .chi_averager import average_chi_files, average_chi_from_database
from .build_chi_database import build_database
from .chi_database_query import ChiDatabaseQuery, PathQuery

__all__ = [
    "load_config", 
    "ConfigError", 
    "process_trajectory", 
    "average_chi_files",
    "average_chi_from_database",
    "build_database",
    "ChiDatabaseQuery",
    "PathQuery"
]
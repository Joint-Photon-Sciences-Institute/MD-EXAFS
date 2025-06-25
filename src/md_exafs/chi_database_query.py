"""
Query functions for the chi(k) database.

This module provides functions to query the chi(k) database for paths
matching specific criteria and retrieve pre-computed chi(k) data.
"""

import sqlite3
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
import numpy as np
import json
from dataclasses import dataclass

from .db_schema import bytes_to_numpy

logger = logging.getLogger(__name__)


@dataclass
class PathQuery:
    """Parameters for querying paths from the database."""
    path_types: Optional[List[str]] = None
    min_reff: Optional[float] = None
    max_reff: Optional[float] = None
    frames: Optional[List[int]] = None
    atom_ids: Optional[List[int]] = None
    nleg: Optional[int] = None


class ChiDatabaseQuery:
    """Query interface for the chi(k) database."""
    
    def __init__(self, db_path: Path):
        """
        Initialize database query interface.
        
        Args:
            db_path: Path to the SQLite database
        """
        if not db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")
            
        self.db_path = db_path
        self._conn = None
    
    def __enter__(self):
        """Context manager entry."""
        self._conn = sqlite3.connect(self.db_path)
        self._conn.row_factory = sqlite3.Row  # Enable column access by name
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self._conn:
            self._conn.close()
    
    def query_paths(self, query: PathQuery) -> List[Dict[str, Any]]:
        """
        Query paths matching the specified criteria.
        
        Args:
            query: PathQuery object with selection criteria
            
        Returns:
            List of dictionaries with path information
        """
        if not self._conn:
            raise RuntimeError("Database connection not established. Use with context manager.")
        
        # Build SQL query
        sql = """
            SELECT p.id, p.frame, p.atom_id, p.path_num, p.nleg, 
                   p.deg, p.reff, p.path_type, p.atom_sequence
            FROM paths p
            WHERE 1=1
        """
        params = []
        
        # Add filters
        if query.path_types:
            placeholders = ','.join('?' * len(query.path_types))
            sql += f" AND p.path_type IN ({placeholders})"
            params.extend(query.path_types)
        
        if query.min_reff is not None:
            sql += " AND p.reff >= ?"
            params.append(query.min_reff)
        
        if query.max_reff is not None:
            sql += " AND p.reff <= ?"
            params.append(query.max_reff)
        
        if query.frames:
            placeholders = ','.join('?' * len(query.frames))
            sql += f" AND p.frame IN ({placeholders})"
            params.extend(query.frames)
        
        if query.atom_ids:
            placeholders = ','.join('?' * len(query.atom_ids))
            sql += f" AND p.atom_id IN ({placeholders})"
            params.extend(query.atom_ids)
        
        if query.nleg is not None:
            sql += " AND p.nleg = ?"
            params.append(query.nleg)
        
        sql += " ORDER BY p.frame, p.atom_id, p.path_num"
        
        # Execute query
        cursor = self._conn.execute(sql, params)
        
        # Convert results to dictionaries
        results = []
        for row in cursor:
            result = dict(row)
            # Parse JSON atom sequence
            result['atom_sequence'] = json.loads(result['atom_sequence'])
            results.append(result)
        
        return results
    
    def get_chi_data(self, path_ids: List[int]) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
        """
        Retrieve chi(k) data for specified path IDs.
        
        Args:
            path_ids: List of path IDs to retrieve
            
        Returns:
            Dictionary mapping path_id to (k_grid, chi_values) tuples
        """
        if not self._conn:
            raise RuntimeError("Database connection not established. Use with context manager.")
        
        if not path_ids:
            return {}
        
        # Query chi data
        placeholders = ','.join('?' * len(path_ids))
        sql = f"""
            SELECT path_id, k_grid, chi_values 
            FROM chi_data 
            WHERE path_id IN ({placeholders})
        """
        
        cursor = self._conn.execute(sql, path_ids)
        
        # Convert results
        chi_data = {}
        for row in cursor:
            path_id = row['path_id']
            k_grid = bytes_to_numpy(row['k_grid'])
            chi_values = bytes_to_numpy(row['chi_values'])
            chi_data[path_id] = (k_grid, chi_values)
        
        return chi_data
    
    def get_paths_with_chi(self, query: PathQuery) -> List[Dict[str, Any]]:
        """
        Query paths and retrieve their chi(k) data in one operation.
        
        Args:
            query: PathQuery object with selection criteria
            
        Returns:
            List of dictionaries with path info and chi data
        """
        # First query paths
        paths = self.query_paths(query)
        
        if not paths:
            return []
        
        # Get path IDs
        path_ids = [p['id'] for p in paths]
        
        # Get chi data
        chi_data = self.get_chi_data(path_ids)
        
        # Combine results
        for path in paths:
            if path['id'] in chi_data:
                k_grid, chi_values = chi_data[path['id']]
                path['k_grid'] = k_grid
                path['chi_values'] = chi_values
            else:
                path['k_grid'] = None
                path['chi_values'] = None
        
        return paths
    
    def get_unique_path_types(self) -> List[str]:
        """
        Get all unique path types in the database.
        
        Returns:
            List of unique path type strings
        """
        if not self._conn:
            raise RuntimeError("Database connection not established. Use with context manager.")
        
        cursor = self._conn.execute("""
            SELECT DISTINCT path_type 
            FROM paths 
            ORDER BY path_type
        """)
        
        return [row[0] for row in cursor]
    
    def get_reff_range(self) -> Tuple[float, float]:
        """
        Get the range of reff values in the database.
        
        Returns:
            Tuple of (min_reff, max_reff)
        """
        if not self._conn:
            raise RuntimeError("Database connection not established. Use with context manager.")
        
        cursor = self._conn.execute("""
            SELECT MIN(reff), MAX(reff) 
            FROM paths
        """)
        
        row = cursor.fetchone()
        return (row[0] or 0.0, row[1] or 0.0)
    
    def get_frame_atom_counts(self) -> Dict[int, int]:
        """
        Get the number of atoms per frame.
        
        Returns:
            Dictionary mapping frame number to atom count
        """
        if not self._conn:
            raise RuntimeError("Database connection not established. Use with context manager.")
        
        cursor = self._conn.execute("""
            SELECT frame, COUNT(DISTINCT atom_id) as atom_count
            FROM atoms
            GROUP BY frame
            ORDER BY frame
        """)
        
        return {row['frame']: row['atom_count'] for row in cursor}
    
    def average_chi_data(self, path_ids: List[int], 
                        k_min: float = 0.0, k_max: float = 20.0, 
                        k_step: float = 0.05) -> Optional[np.ndarray]:
        """
        Average chi(k) data for specified paths and interpolate to standard grid.
        
        Args:
            path_ids: List of path IDs to average
            k_min: Minimum k value for output grid
            k_max: Maximum k value for output grid
            k_step: Step size for output grid
            
        Returns:
            Numpy array with columns [k, chi_average] or None if no data
        """
        if not path_ids:
            return None
        
        # Get chi data
        chi_data = self.get_chi_data(path_ids)
        
        if not chi_data:
            return None
        
        # Create standard k-grid
        standard_k = np.arange(k_min, k_max + k_step, k_step)
        
        # Interpolate all chi data to standard grid
        interpolated_chi = []
        
        for path_id, (k_grid, chi_values) in chi_data.items():
            # Interpolate to standard grid
            chi_interp = np.interp(standard_k, k_grid, chi_values, left=0, right=0)
            interpolated_chi.append(chi_interp)
        
        # Average the interpolated chi values
        if interpolated_chi:
            chi_average = np.mean(interpolated_chi, axis=0)
            return np.column_stack((standard_k, chi_average))
        else:
            return None
    
    def sum_chi_data(self, path_ids: List[int],
                     k_min: float = 0.0, k_max: float = 20.0,
                     k_step: float = 0.05) -> Optional[np.ndarray]:
        """
        Sum chi(k) data for specified paths and interpolate to standard grid.
        
        Args:
            path_ids: List of path IDs to sum
            k_min: Minimum k value for output grid
            k_max: Maximum k value for output grid
            k_step: Step size for output grid
            
        Returns:
            Numpy array with columns [k, chi_sum] or None if no data
        """
        if not path_ids:
            return None
        
        # Get chi data
        chi_data = self.get_chi_data(path_ids)
        
        if not chi_data:
            return None
        
        # Create standard k-grid
        standard_k = np.arange(k_min, k_max + k_step, k_step)
        
        # Interpolate all chi data to standard grid and sum
        chi_sum = np.zeros_like(standard_k)
        
        for path_id, (k_grid, chi_values) in chi_data.items():
            # Interpolate to standard grid
            chi_interp = np.interp(standard_k, k_grid, chi_values, left=0, right=0)
            chi_sum += chi_interp
        
        return np.column_stack((standard_k, chi_sum))


def query_and_average_multipath(db_path: Path, 
                               path_types: List[str],
                               max_distances: List[Optional[float]],
                               frames: Optional[List[int]] = None) -> Optional[np.ndarray]:
    """
    Query database for multipath selection and return averaged chi(k).
    
    Args:
        db_path: Path to the database
        path_types: List of path types to include
        max_distances: List of maximum distances for each path type
        frames: Optional list of frames to include
        
    Returns:
        Averaged chi(k) data as numpy array with columns [k, chi]
    """
    with ChiDatabaseQuery(db_path) as db:
        all_path_ids = []
        
        # Query each path type
        for path_type, max_dist in zip(path_types, max_distances):
            query = PathQuery(
                path_types=[path_type],
                max_reff=max_dist,
                frames=frames
            )
            
            paths = db.query_paths(query)
            path_ids = [p['id'] for p in paths]
            all_path_ids.extend(path_ids)
            
            logger.info(f"Found {len(path_ids)} paths of type '{path_type}' "
                       f"with reff <= {max_dist}")
        
        if not all_path_ids:
            logger.warning("No paths found matching criteria")
            return None
        
        # Average chi data
        logger.info(f"Averaging {len(all_path_ids)} total paths")
        return db.average_chi_data(all_path_ids)
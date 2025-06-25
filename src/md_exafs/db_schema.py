"""
Database schema and helper functions for MD-EXAFS chi(k) data caching.

This module provides SQLite database functionality for storing pre-computed
FEFF path metadata and chi(k) data to accelerate multipath averaging operations.
"""

import sqlite3
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import json
import logging

logger = logging.getLogger(__name__)


def create_database(db_path: Path) -> None:
    """
    Create the SQLite database with the required schema.
    
    Args:
        db_path: Path to the database file to create
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Enable foreign keys
    cursor.execute("PRAGMA foreign_keys = ON")
    
    # Create tables
    cursor.executescript("""
        -- Store path metadata
        CREATE TABLE IF NOT EXISTS paths (
            id INTEGER PRIMARY KEY,
            frame INTEGER NOT NULL,
            atom_id INTEGER NOT NULL,
            path_num INTEGER NOT NULL,
            nleg INTEGER NOT NULL,
            deg REAL NOT NULL,
            reff REAL NOT NULL,
            path_type TEXT NOT NULL,
            atom_sequence TEXT NOT NULL,  -- JSON array of atoms
            UNIQUE(frame, atom_id, path_num)
        );
        
        -- Store pre-computed chi data
        CREATE TABLE IF NOT EXISTS chi_data (
            id INTEGER PRIMARY KEY,
            path_id INTEGER NOT NULL,
            k_grid BLOB NOT NULL,  -- Numpy array as bytes
            chi_values BLOB NOT NULL,  -- Numpy array as bytes
            FOREIGN KEY(path_id) REFERENCES paths(id) ON DELETE CASCADE
        );
        
        -- Store atom information
        CREATE TABLE IF NOT EXISTS atoms (
            id INTEGER PRIMARY KEY,
            frame INTEGER NOT NULL,
            atom_id INTEGER NOT NULL,
            element TEXT NOT NULL,
            position_x REAL,
            position_y REAL,
            position_z REAL,
            UNIQUE(frame, atom_id)
        );
        
        -- Create indexes for fast querying
        CREATE INDEX IF NOT EXISTS idx_path_type ON paths(path_type);
        CREATE INDEX IF NOT EXISTS idx_reff ON paths(reff);
        CREATE INDEX IF NOT EXISTS idx_frame_atom ON paths(frame, atom_id);
        CREATE INDEX IF NOT EXISTS idx_chi_path ON chi_data(path_id);
        
        -- Metadata table for database versioning
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
    """)
    
    # Store database version
    cursor.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
                   ("schema_version", "1.0"))
    cursor.execute("INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
                   ("created_timestamp", "datetime('now')"))
    
    conn.commit()
    conn.close()
    logger.info(f"Created database at {db_path}")


def validate_database(db_path: Path) -> bool:
    """
    Validate that the database exists and has the correct schema.
    
    Args:
        db_path: Path to the database file
        
    Returns:
        True if database is valid, False otherwise
    """
    if not db_path.exists():
        return False
        
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check for required tables
        required_tables = ['paths', 'chi_data', 'atoms', 'metadata']
        cursor.execute("""
            SELECT name FROM sqlite_master 
            WHERE type='table' AND name IN ({})
        """.format(','.join('?' * len(required_tables))), required_tables)
        
        found_tables = {row[0] for row in cursor.fetchall()}
        
        conn.close()
        
        return len(found_tables) == len(required_tables)
        
    except sqlite3.Error as e:
        logger.error(f"Database validation error: {e}")
        return False


def numpy_to_bytes(arr: np.ndarray) -> bytes:
    """
    Convert numpy array to bytes for database storage.
    
    Args:
        arr: Numpy array to convert
        
    Returns:
        Bytes representation of the array
    """
    return arr.tobytes()


def bytes_to_numpy(data: bytes, dtype: np.dtype = np.float64) -> np.ndarray:
    """
    Convert bytes back to numpy array.
    
    Args:
        data: Bytes data from database
        dtype: Data type of the array
        
    Returns:
        Reconstructed numpy array
    """
    return np.frombuffer(data, dtype=dtype)


def insert_path_metadata(conn: sqlite3.Connection, 
                        frame: int, 
                        atom_id: int, 
                        path_num: int,
                        path_info: Dict[str, Any]) -> int:
    """
    Insert path metadata into the database.
    
    Args:
        conn: Database connection
        frame: Frame number
        atom_id: Atom ID
        path_num: Path number
        path_info: Dictionary containing path metadata
        
    Returns:
        ID of the inserted path record
    """
    cursor = conn.cursor()
    
    # Convert atom sequence to JSON
    atom_sequence_json = json.dumps(path_info['atoms'])
    
    cursor.execute("""
        INSERT OR REPLACE INTO paths 
        (frame, atom_id, path_num, nleg, deg, reff, path_type, atom_sequence)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, (frame, atom_id, path_num, 
          path_info['nleg'], path_info['deg'], path_info['reff'],
          path_info['path_type'], atom_sequence_json))
    
    return cursor.lastrowid


def insert_chi_data(conn: sqlite3.Connection,
                   path_id: int,
                   k_grid: np.ndarray,
                   chi_values: np.ndarray) -> None:
    """
    Insert pre-computed chi(k) data into the database.
    
    Args:
        conn: Database connection
        path_id: ID of the associated path record
        k_grid: k-space grid points
        chi_values: chi(k) values
    """
    cursor = conn.cursor()
    
    k_bytes = numpy_to_bytes(k_grid)
    chi_bytes = numpy_to_bytes(chi_values)
    
    cursor.execute("""
        INSERT OR REPLACE INTO chi_data (path_id, k_grid, chi_values)
        VALUES (?, ?, ?)
    """, (path_id, k_bytes, chi_bytes))


def insert_atom_info(conn: sqlite3.Connection,
                    frame: int,
                    atom_id: int,
                    element: str,
                    position: Optional[Tuple[float, float, float]] = None) -> None:
    """
    Insert atom information into the database.
    
    Args:
        conn: Database connection
        frame: Frame number
        atom_id: Atom ID
        element: Element symbol
        position: Optional (x, y, z) coordinates
    """
    cursor = conn.cursor()
    
    if position:
        x, y, z = position
    else:
        x = y = z = None
        
    cursor.execute("""
        INSERT OR REPLACE INTO atoms 
        (frame, atom_id, element, position_x, position_y, position_z)
        VALUES (?, ?, ?, ?, ?, ?)
    """, (frame, atom_id, element, x, y, z))


def get_database_stats(db_path: Path) -> Dict[str, Any]:
    """
    Get statistics about the database contents.
    
    Args:
        db_path: Path to the database file
        
    Returns:
        Dictionary with database statistics
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    stats = {}
    
    # Get counts
    cursor.execute("SELECT COUNT(*) FROM paths")
    stats['total_paths'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM chi_data")
    stats['total_chi_data'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT frame) FROM paths")
    stats['total_frames'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(DISTINCT atom_id) FROM atoms")
    stats['total_atoms'] = cursor.fetchone()[0]
    
    # Get path type distribution
    cursor.execute("""
        SELECT path_type, COUNT(*) 
        FROM paths 
        GROUP BY path_type
        ORDER BY COUNT(*) DESC
    """)
    stats['path_types'] = dict(cursor.fetchall())
    
    # Get database size
    cursor.execute("SELECT page_count * page_size FROM pragma_page_count(), pragma_page_size()")
    stats['size_bytes'] = cursor.fetchone()[0]
    
    # Get metadata
    cursor.execute("SELECT key, value FROM metadata")
    stats['metadata'] = dict(cursor.fetchall())
    
    conn.close()
    
    return stats


def clear_database(db_path: Path) -> None:
    """
    Clear all data from the database while preserving schema.
    
    Args:
        db_path: Path to the database file
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Delete all data
    cursor.execute("DELETE FROM chi_data")
    cursor.execute("DELETE FROM paths")
    cursor.execute("DELETE FROM atoms")
    
    # Update metadata
    cursor.execute("UPDATE metadata SET value = datetime('now') WHERE key = 'last_cleared'")
    
    # Vacuum to reclaim space
    cursor.execute("VACUUM")
    
    conn.commit()
    conn.close()
    
    logger.info(f"Cleared all data from database at {db_path}")


def optimize_database(db_path: Path) -> None:
    """
    Optimize database for better query performance.
    
    Args:
        db_path: Path to the database file
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Analyze tables for query optimization
    cursor.execute("ANALYZE")
    
    # Vacuum to defragment
    cursor.execute("VACUUM")
    
    conn.commit()
    conn.close()
    
    logger.info(f"Optimized database at {db_path}")
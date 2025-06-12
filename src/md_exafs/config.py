"""Configuration parser and validator for MD-EXAFS."""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


class ConfigError(Exception):
    """Raised when configuration is invalid or missing required fields."""
    pass


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load and validate configuration from TOML file.
    
    Args:
        config_path: Path to TOML configuration file
        
    Returns:
        Validated configuration dictionary
        
    Raises:
        ConfigError: If configuration is invalid or missing required fields
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise ConfigError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_path, "rb") as f:
            config = tomllib.load(f)
    except Exception as e:
        raise ConfigError(f"Failed to parse TOML file: {e}")
    
    # Validate configuration
    _validate_config(config)
    
    # Set default FEFF executable path if not provided
    if "feff_executable" not in config.get("feff", {}):
        # Use bundled executable by default
        bundled_feff = Path(__file__).parent / "feff8_lite" / "feff85L.exe"
        config.setdefault("feff", {})["feff_executable"] = str(bundled_feff)
    
    return config


def _validate_config(config: Dict[str, Any]) -> None:
    """Validate configuration has all required fields."""
    # Determine config type based on presence of sections
    if "averaging" in config:
        _validate_averaging_config(config)
    else:
        _validate_processing_config(config)


def _validate_processing_config(config: Dict[str, Any]) -> None:
    """Validate processing configuration."""
    # Check required top-level sections
    required_sections = ["system", "lattice", "atoms", "processing", "feff"]
    for section in required_sections:
        if section not in config:
            raise ConfigError(f"Missing required section: [{section}]")
    
    # Validate system section
    system = config["system"]
    required_system_fields = ["name", "trajectory_file", "output_directory"]
    for field in required_system_fields:
        if field not in system:
            raise ConfigError(f"Missing required field in [system]: {field}")
    
    # Check trajectory file exists
    trajectory_path = Path(system["trajectory_file"])
    if not trajectory_path.exists():
        raise ConfigError(f"Trajectory file not found: {trajectory_path}")
    
    # Validate lattice section
    lattice = config["lattice"]
    required_lattice_fields = ["a", "b", "c", "pbc"]
    for field in required_lattice_fields:
        if field not in lattice:
            raise ConfigError(f"Missing required field in [lattice]: {field}")
    
    # Validate lattice vectors are 3D
    for vec_name in ["a", "b", "c"]:
        vec = lattice[vec_name]
        if not isinstance(vec, list) or len(vec) != 3:
            raise ConfigError(f"Lattice vector {vec_name} must be a list of 3 numbers")
        if not all(isinstance(x, (int, float)) for x in vec):
            raise ConfigError(f"Lattice vector {vec_name} must contain only numbers")
    
    # Validate PBC
    pbc = lattice["pbc"]
    if not isinstance(pbc, list) or len(pbc) != 3:
        raise ConfigError("pbc must be a list of 3 boolean values")
    if not all(isinstance(x, bool) for x in pbc):
        raise ConfigError("pbc must contain only boolean values")
    
    # Validate atoms section
    atoms = config["atoms"]
    if not atoms:
        raise ConfigError("[atoms] section cannot be empty")
    
    # Validate processing section
    processing = config["processing"]
    required_processing_fields = [
        "start_frame", "end_frame", "frame_step", 
        "atoms_per_frame", "cutoff_radius"
    ]
    for field in required_processing_fields:
        if field not in processing:
            raise ConfigError(f"Missing required field in [processing]: {field}")
    
    # Validate frame range
    if processing["start_frame"] >= processing["end_frame"]:
        raise ConfigError("start_frame must be less than end_frame")
    
    if processing["frame_step"] <= 0:
        raise ConfigError("frame_step must be positive")
    
    if processing["atoms_per_frame"] <= 0:
        raise ConfigError("atoms_per_frame must be positive")
    
    if processing["cutoff_radius"] <= 0:
        raise ConfigError("cutoff_radius must be positive")
    
    # Validate FEFF section
    feff = config["feff"]
    required_feff_fields = ["header", "absorbing_element"]
    for field in required_feff_fields:
        if field not in feff:
            raise ConfigError(f"Missing required field in [feff]: {field}")
    
    # Validate absorbing element is in atoms
    absorbing_element = feff["absorbing_element"]
    if absorbing_element not in atoms:
        raise ConfigError(f"Absorbing element '{absorbing_element}' not found in [atoms] section")


def _validate_averaging_config(config: Dict[str, Any]) -> None:
    """Validate averaging configuration."""
    if "averaging" not in config:
        raise ConfigError("Missing required section: [averaging]")
    
    averaging = config["averaging"]
    required_fields = ["input_directory", "frame_range", "output_file"]
    
    for field in required_fields:
        if field not in averaging:
            raise ConfigError(f"Missing required field in [averaging]: {field}")
    
    # Validate frame range
    frame_range = averaging["frame_range"]
    if not isinstance(frame_range, list) or len(frame_range) != 2:
        raise ConfigError("frame_range must be [start, end]")
    if frame_range[0] >= frame_range[1]:
        raise ConfigError("frame_range: start must be less than end")
    
    # Check input directory exists
    input_dir = Path(averaging["input_directory"])
    if not input_dir.exists():
        raise ConfigError(f"Input directory not found: {input_dir}")


def get_element_data() -> Dict[int, Tuple[str, str]]:
    """
    Get periodic table data for element support.
    
    Returns:
        Dictionary mapping atomic number to (symbol, name)
    """
    return {
        1: ("H", "Hydrogen"), 2: ("He", "Helium"),
        3: ("Li", "Lithium"), 4: ("Be", "Beryllium"), 5: ("B", "Boron"),
        6: ("C", "Carbon"), 7: ("N", "Nitrogen"), 8: ("O", "Oxygen"),
        9: ("F", "Fluorine"), 10: ("Ne", "Neon"),
        11: ("Na", "Sodium"), 12: ("Mg", "Magnesium"), 13: ("Al", "Aluminum"),
        14: ("Si", "Silicon"), 15: ("P", "Phosphorus"), 16: ("S", "Sulfur"),
        17: ("Cl", "Chlorine"), 18: ("Ar", "Argon"),
        19: ("K", "Potassium"), 20: ("Ca", "Calcium"), 21: ("Sc", "Scandium"),
        22: ("Ti", "Titanium"), 23: ("V", "Vanadium"), 24: ("Cr", "Chromium"),
        25: ("Mn", "Manganese"), 26: ("Fe", "Iron"), 27: ("Co", "Cobalt"),
        28: ("Ni", "Nickel"), 29: ("Cu", "Copper"), 30: ("Zn", "Zinc"),
        31: ("Ga", "Gallium"), 32: ("Ge", "Germanium"), 33: ("As", "Arsenic"),
        34: ("Se", "Selenium"), 35: ("Br", "Bromine"), 36: ("Kr", "Krypton"),
        37: ("Rb", "Rubidium"), 38: ("Sr", "Strontium"), 39: ("Y", "Yttrium"),
        40: ("Zr", "Zirconium"), 41: ("Nb", "Niobium"), 42: ("Mo", "Molybdenum"),
        43: ("Tc", "Technetium"), 44: ("Ru", "Ruthenium"), 45: ("Rh", "Rhodium"),
        46: ("Pd", "Palladium"), 47: ("Ag", "Silver"), 48: ("Cd", "Cadmium"),
        49: ("In", "Indium"), 50: ("Sn", "Tin"), 51: ("Sb", "Antimony"),
        52: ("Te", "Tellurium"), 53: ("I", "Iodine"), 54: ("Xe", "Xenon"),
        55: ("Cs", "Cesium"), 56: ("Ba", "Barium"), 57: ("La", "Lanthanum"),
        58: ("Ce", "Cerium"), 59: ("Pr", "Praseodymium"), 60: ("Nd", "Neodymium"),
        61: ("Pm", "Promethium"), 62: ("Sm", "Samarium"), 63: ("Eu", "Europium"),
        64: ("Gd", "Gadolinium"), 65: ("Tb", "Terbium"), 66: ("Dy", "Dysprosium"),
        67: ("Ho", "Holmium"), 68: ("Er", "Erbium"), 69: ("Tm", "Thulium"),
        70: ("Yb", "Ytterbium"), 71: ("Lu", "Lutetium"), 72: ("Hf", "Hafnium"),
        73: ("Ta", "Tantalum"), 74: ("W", "Tungsten"), 75: ("Re", "Rhenium"),
        76: ("Os", "Osmium"), 77: ("Ir", "Iridium"), 78: ("Pt", "Platinum"),
        79: ("Au", "Gold"), 80: ("Hg", "Mercury"), 81: ("Tl", "Thallium"),
        82: ("Pb", "Lead"), 83: ("Bi", "Bismuth"), 84: ("Po", "Polonium"),
        85: ("At", "Astatine"), 86: ("Rn", "Radon"), 87: ("Fr", "Francium"),
        88: ("Ra", "Radium"), 89: ("Ac", "Actinium"), 90: ("Th", "Thorium"),
        91: ("Pa", "Protactinium"), 92: ("U", "Uranium"), 93: ("Np", "Neptunium"),
        94: ("Pu", "Plutonium"), 95: ("Am", "Americium"), 96: ("Cm", "Curium"),
        97: ("Bk", "Berkelium"), 98: ("Cf", "Californium"), 99: ("Es", "Einsteinium"),
        100: ("Fm", "Fermium"), 101: ("Md", "Mendelevium"), 102: ("No", "Nobelium"),
        103: ("Lr", "Lawrencium"),
    }


def get_atomic_number(symbol: str) -> int:
    """
    Get atomic number from element symbol.
    
    Args:
        symbol: Element symbol (e.g., "U", "O")
        
    Returns:
        Atomic number
        
    Raises:
        ConfigError: If symbol is not recognized
    """
    element_data = get_element_data()
    for atomic_num, (sym, _) in element_data.items():
        if sym == symbol:
            return atomic_num
    
    raise ConfigError(f"Unknown element symbol: {symbol}")
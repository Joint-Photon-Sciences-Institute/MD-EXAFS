"""Check MD convergence by plotting energy data from .ener files."""

import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def parse_ener_file(filepath: Path) -> Dict[str, np.ndarray]:
    """Parse CP2K .ener file and extract energy data.
    
    Args:
        filepath: Path to .ener file
        
    Returns:
        Dictionary with arrays for each column
    """
    data = {
        'step': [],
        'time': [],
        'kinetic': [],
        'temperature': [],
        'potential': [],
        'conserved': [],
        'used_time': []
    }
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 7:
                data['step'].append(int(parts[0]))
                data['time'].append(float(parts[1]))
                data['kinetic'].append(float(parts[2]))
                data['temperature'].append(float(parts[3]))
                data['potential'].append(float(parts[4]))
                data['conserved'].append(float(parts[5]))
                data['used_time'].append(float(parts[6]))
    
    # Convert to numpy arrays
    for key in data:
        data[key] = np.array(data[key])
    
    return data


def calculate_running_average(data: np.ndarray, window_size: int = 100) -> np.ndarray:
    """Calculate running average with given window size."""
    if len(data) < window_size:
        window_size = max(1, len(data) // 10)
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


def plot_convergence(data: Dict[str, np.ndarray], output_path: Path) -> None:
    """Create convergence plots for MD simulation.
    
    Args:
        data: Dictionary with energy data arrays
        output_path: Path to save the plot
    """
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)
    
    # Temperature plot
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(data['time'], data['temperature'], alpha=0.6, label='Instantaneous')
    temp_avg = calculate_running_average(data['temperature'])
    time_avg = data['time'][:len(temp_avg)]
    ax1.plot(time_avg, temp_avg, 'r-', linewidth=2, label='Running average')
    ax1.set_xlabel('Time (fs)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Temperature Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Total energy (conserved quantity) plot
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(data['time'], data['conserved'], alpha=0.6)
    ax2.set_xlabel('Time (fs)')
    ax2.set_ylabel('Conserved Energy (a.u.)')
    ax2.set_title('Total Energy Conservation')
    ax2.grid(True, alpha=0.3)
    
    # Kinetic energy plot
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(data['time'], data['kinetic'], alpha=0.6, label='Instantaneous')
    kin_avg = calculate_running_average(data['kinetic'])
    ax3.plot(time_avg, kin_avg, 'r-', linewidth=2, label='Running average')
    ax3.set_xlabel('Time (fs)')
    ax3.set_ylabel('Kinetic Energy (a.u.)')
    ax3.set_title('Kinetic Energy Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Potential energy plot
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(data['time'], data['potential'], alpha=0.6, label='Instantaneous')
    pot_avg = calculate_running_average(data['potential'])
    ax4.plot(time_avg, pot_avg, 'r-', linewidth=2, label='Running average')
    ax4.set_xlabel('Time (fs)')
    ax4.set_ylabel('Potential Energy (a.u.)')
    ax4.set_title('Potential Energy Evolution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Energy drift plot (conserved energy - initial value)
    ax5 = fig.add_subplot(gs[2, :])
    energy_drift = (data['conserved'] - data['conserved'][0]) * 627.509  # Convert to kcal/mol
    ax5.plot(data['time'], energy_drift)
    ax5.set_xlabel('Time (fs)')
    ax5.set_ylabel('Energy Drift (kcal/mol)')
    ax5.set_title('Total Energy Drift from Initial Value')
    ax5.grid(True, alpha=0.3)
    
    # Add overall title
    fig.suptitle('MD Simulation Convergence Analysis', fontsize=14, fontweight='bold')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Print convergence statistics
    print(f"\nConvergence Statistics:")
    print(f"  Temperature: {np.mean(data['temperature']):.2f} ± {np.std(data['temperature']):.2f} K")
    print(f"  Energy drift: {energy_drift[-1]:.4f} kcal/mol")
    print(f"  Relative energy drift: {(energy_drift[-1] / (data['conserved'][0] * 627.509)) * 100:.4e} %")
    print(f"  Simulation time: {data['time'][-1]:.1f} fs")
    print(f"\nPlot saved to: {output_path}")


def main():
    """Main function for convergence checking CLI."""
    parser = argparse.ArgumentParser(
        description='Check MD simulation convergence by plotting energy data from .ener files'
    )
    parser.add_argument('--check_convergence', action='store_true',
                        help='Enable convergence checking mode')
    parser.add_argument('--input', type=str, required=True,
                        help='Path to .ener file')
    parser.add_argument('--output', type=str, required=True,
                        help='Path for output plot (PNG)')
    
    args = parser.parse_args()
    
    if not args.check_convergence:
        parser.error("--check_convergence flag is required")
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    # Parse energy file
    print(f"Reading energy file: {input_path}")
    data = parse_ener_file(input_path)
    
    if len(data['step']) == 0:
        raise ValueError("No data found in energy file")
    
    # Create convergence plots
    print("Creating convergence plots...")
    plot_convergence(data, output_path)


if __name__ == '__main__':
    main()
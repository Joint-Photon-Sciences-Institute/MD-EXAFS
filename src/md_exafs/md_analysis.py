"""MD analysis tool for convergence checking and RDF calculations."""

import argparse
import sys
from pathlib import Path

# Import the specific analysis modules
from .convergence_check import parse_ener_file, plot_convergence
from .rdf_analysis import main as rdf_main


def convergence_analysis(args):
    """Handle convergence checking."""
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


def main():
    """Main entry point for MD analysis tools."""
    parser = argparse.ArgumentParser(
        description='MD analysis tools for convergence checking and RDF calculations',
        epilog='Use either --check_convergence or --rdf to specify the analysis type'
    )
    
    # Convergence checking arguments
    parser.add_argument('--check_convergence', action='store_true',
                       help='Enable convergence checking mode')
    parser.add_argument('--input', type=str,
                       help='Path to .ener file (for convergence checking)')
    parser.add_argument('--output', type=str,
                       help='Path for output plot (PNG) (for convergence checking)')
    
    # RDF analysis arguments
    parser.add_argument('--rdf', type=str,
                       help='Path to RDF configuration TOML file')
    
    args = parser.parse_args()
    
    # Check that exactly one mode is selected
    if args.check_convergence and args.rdf:
        parser.error("Please specify either --check_convergence or --rdf, not both")
    
    if not args.check_convergence and not args.rdf:
        parser.error("Please specify either --check_convergence or --rdf")
    
    # Handle convergence checking
    if args.check_convergence:
        if not args.input or not args.output:
            parser.error("--check_convergence requires --input and --output arguments")
        convergence_analysis(args)
    
    # Handle RDF analysis
    elif args.rdf:
        # Pass the argument to the RDF analysis module
        sys.argv = ['md-exafs-md', '--rdf', args.rdf]
        rdf_main()


if __name__ == '__main__':
    main()
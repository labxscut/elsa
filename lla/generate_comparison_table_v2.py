'''
No parallel processing in this version.
Only consider D=0 in this version.
'''

#!/usr/bin/env python

import numpy as np
import pandas as pd
import tempfile
import sys
import argparse
import time
import os
import traceback
import logging

# Ensure consistent imports
try:
    from lsa import lsalib
    from lla import lla_compute, lla_sim
    print("Using package imports for lsalib, lla_compute, and lla_sim", file=sys.stderr)
except ImportError:
    try:
        import lsalib
        import lla_compute
        import lla_sim
        print("Using direct imports for lsalib, lla_compute, and lla_sim", file=sys.stderr)
    except ImportError:
        print("ERROR: Could not import required modules. Please ensure lsa and lla packages are installed.", file=sys.stderr)
        sys.exit(1)

def setup_logging(debug=False):
    """Setup logging configuration"""
    level = logging.DEBUG if debug else logging.INFO
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"comparison_table_v2_{timestamp}.log")
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging initialized. Log file: {log_file}")

def generate_simulated_data(n_points, n_series=100, rep_num=1, debug=False):
    """Generate simulated time series data for different sample sizes
    
    Args:
        n_points: Number of time points
        n_series: Number of series to generate
        rep_num: Number of replicates per time point
        debug: Whether to print debug information
        
    Returns:
        Path to the generated data file
    """
    if debug:
        logging.debug(f"Generating simulated data with n_points={n_points}, n_series={n_series}, rep_num={rep_num}")
    
    # Create a named temporary file that won't be deleted when closed
    output_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
    output_path = output_file.name
    output_file.close()  # Close the file so lla_sim can write to it
    
    # Save original argv
    original_argv = sys.argv
    
    # Set up arguments for lla_sim
    sys.argv = [
        'lla_sim',
        '-L', str(n_points),
        '-N', str(n_series),
        '-R', str(rep_num),
        '-M', 'idn,0,1',
        '-O', output_path
    ]
    
    try:
        if debug:
            logging.debug(f"Running lla_sim with args: {' '.join(sys.argv)}")
        
        # Call lla_sim.main()
        lla_sim.main()
        
        # Verify the file was created and has content
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            raise FileNotFoundError(f"Failed to generate data file or file is empty: {output_path}")
            
        if debug:
            logging.debug(f"Successfully generated data file: {output_path}")
            
        return output_path
    except Exception as e:
        logging.error(f"Error in generate_simulated_data: {str(e)}")
        if debug:
            logging.debug("Traceback:", exc_info=True)
        if os.path.exists(output_path):
            try:
                os.unlink(output_path)
            except:
                pass
        raise
    finally:
        # Restore original argv
        sys.argv = original_argv

def calculate_theoretical_pvalues(x_values, Rmax=None, Dmax=0, precision=0.001, debug=False):
    """Calculate theoretical P-values with adaptive parameters"""
    if Rmax is None:
        Rmax = max(x_values) * 2
    Rmax = min(max(Rmax, 10), 50)  # Ensure Rmax is within bounds
    
    if debug:
        logging.debug(f"Calculating theoretical P-values with Rmax={Rmax}, Dmax={Dmax}")
    
    P_table = lsalib.theoPvalue(Rmax=Rmax, Dmax=Dmax, precision=precision)
    return [lsalib.readPvalue(P_table, x, 1) for x in x_values]

def calculate_simulated_pvalues(data_file, n_points, x_value, delay=0, precision=1000, debug=False):
    """Calculate simulated P-values for a single x value"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as output_file:
        output_path = output_file.name
    
    if debug:
        logging.debug(f"Calculating simulated P-value for x={x_value}")
    
    try:
        sys.argv = [
            'lla_compute',
            data_file,
            output_path,
            "-d", str(delay),
            "-p", "perm",
            "-x", str(precision),
            "-r", "1",
            "-s", str(n_points),
            "-m", "0.5",
            "-t", "simple",
            "-f", "linear",
            "-n", "pnz"
        ]
        
        if debug:
            logging.debug(f"Running lla_compute with args: {' '.join(sys.argv)}")
        
        lla_compute.main()
        
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            logging.error("Output file not found or empty")
            return 0.0
            
        results = pd.read_csv(output_path, sep='\t', comment='#')
        
        if 'LA' not in results.columns:
            logging.error(f"LA column not found. Available columns: {results.columns.tolist()}")
            return 0.0
            
        la_values = results['LA'].abs()
        filtered_results = results[la_values >= x_value]
        
        if len(results) == 0:
            logging.warning("No results found in the output file")
            return 0.0
            
        p_value = len(filtered_results) / len(results)
        
        if debug:
            logging.debug(f"Total triplets: {len(results)}")
            logging.debug(f"Filtered triplets: {len(filtered_results)}")
            logging.debug(f"LA value range: [{la_values.min():.4f}, {la_values.max():.4f}]")
            logging.debug(f"Threshold x_value: {x_value}")
            logging.debug(f"Calculated P-value: {p_value}")
            
        return p_value
        
    except Exception as e:
        logging.error(f"Error in calculate_simulated_pvalues: {str(e)}")
        if debug:
            logging.debug("Traceback:", exc_info=True)
        return 0.0
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)

def calculate_all_pvalues(data_file, n_points, x_values, debug=False):
    """Calculate P-values for all x values sequentially"""
    sim_pvalues = []
    for x in x_values:
        if debug:
            logging.debug(f"Processing x value: {x}")
        p_value = calculate_simulated_pvalues(data_file, n_points, x, debug=debug)
        sim_pvalues.append(p_value)
    return sim_pvalues

def main():
    """Main function to run the comparison analysis"""
    parser = argparse.ArgumentParser(description="Generate comparison table of theoretical vs simulated P-values")
    parser.add_argument("--xmin", type=float, default=2.0)
    parser.add_argument("--xmax", type=float, default=5.2)
    parser.add_argument("--xstep", type=float, default=0.2)
    parser.add_argument("--sample_sizes", type=int, nargs="+", default=[10, 20, 30])
    parser.add_argument("--n_series", type=int, default=100)
    parser.add_argument("--debug", action="store_true")
    
    args = parser.parse_args()
    setup_logging(args.debug)
    logging.info("Starting comparison analysis")
    
    try:
        x_values = np.arange(args.xmin, args.xmax, args.xstep)
        theo_pvalues = calculate_theoretical_pvalues(x_values, debug=args.debug)
        
        results = pd.DataFrame({
            'x': x_values,
            'Theory': theo_pvalues
        })
        
        for n in args.sample_sizes:
            logging.info(f"Processing sample size: {n}")
            try:
                data_file = generate_simulated_data(n, args.n_series, debug=args.debug)
                sim_pvalues = calculate_all_pvalues(data_file, n, x_values, debug=args.debug)
                results[str(n)] = sim_pvalues
                logging.info(f"Completed processing for sample size {n}")
            finally:
                if os.path.exists(data_file):
                    os.unlink(data_file)
        
        output_file = "comparison_table_v2.txt"
        results.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Results saved to {output_file}")
        
    except Exception as e:
        logging.error(f"Error in main: {str(e)}")
        if args.debug:
            logging.debug("Traceback:", exc_info=True)
        raise

if __name__ == "__main__":
    main()
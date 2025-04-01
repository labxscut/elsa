'''
Test the calculate_simulated_pvalues function.
Parameters are set as follows:
n_points = 10
n_series = 200
x_values = [1.0, 1.5, 2.0, 2.5, 3.0]
debug = True

Output files:
1. Logs:
   - logs/test_simulated_pvalues_20250328_*.log (timestamp-based log file)
   - Contains detailed execution information and debug messages

2. Results:
   - results/test_simulated_pvalues_20250328_*.txt (summary of all results)
   - Contains x values, p-values, and computation times

3. Debug files (if errors occur):
   - debug_output_*.txt (temporary files for debugging)
'''

#!/usr/bin/env python

import sys
import os
import logging
import tempfile
import pandas as pd
import numpy as np
import time

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lla.generate_comparison_table_v2 import calculate_simulated_pvalues, generate_simulated_data

def setup_logging():
    """Setup logging configuration with timestamp-based filename"""
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"test_simulated_pvalues_20250328_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging initialized. Log file: {log_file}")
    return timestamp

def test_simulated_pvalues():
    """Test the calculate_simulated_pvalues function with different parameters"""
    timestamp = setup_logging()
    
    # Test parameters
    n_points = 10  # 10 time points
    n_series = 200  # 200 series
    x_values = [1.0, 1.5, 2.0, 2.5, 3.0]  # x values from 1.0 to 3.0 with step 0.5
    
    # Create results directory and file
    results_dir = "results"
    os.makedirs(results_dir, exist_ok=True)
    results_file = os.path.join(results_dir, f"test_simulated_pvalues_20250328_{timestamp}.txt")
    
    # Initialize results list
    results = []
    
    try:
        # Generate test data
        logging.info(f"Generating test data with n_points={n_points}, n_series={n_series}")
        data_file = generate_simulated_data(n_points, n_series, debug=True)
        
        # Test different x values
        for x in x_values:
            logging.info(f"\nTesting with x={x}")
            start_time = time.time()
            
            p_value = calculate_simulated_pvalues(
                data_file=data_file,
                n_points=n_points,
                x_value=x,
                debug=True
            )
            
            end_time = time.time()
            elapsed_time = end_time - start_time
            
            # Store result
            results.append({
                'x_value': x,
                'p_value': p_value,
                'time_elapsed': elapsed_time
            })
            
            # Log result
            logging.info(f"Result for x={x}: p_value={p_value}")
            logging.info(f"Time elapsed: {elapsed_time:.2f} seconds")
        
        # Save results to file
        results_df = pd.DataFrame(results)
        results_df.to_csv(results_file, sep='\t', index=False)
        logging.info(f"\nResults saved to: {results_file}")
            
    except Exception as e:
        logging.error(f"Error in test: {str(e)}")
        logging.error("Traceback:", exc_info=True)
    finally:
        # Clean up
        if 'data_file' in locals() and os.path.exists(data_file):
            os.unlink(data_file)

if __name__ == "__main__":
    test_simulated_pvalues() 
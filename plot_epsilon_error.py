#!/usr/bin/env python3
"""
Plot epsilon error analysis results from CSV file.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def plot_epsilon_error(csv_file, output_file=None, log_scale=True):
    """
    Plot epsilon vs error from CSV file.
    
    Args:
        csv_file: Path to CSV file with epsilon and error columns
        output_file: Optional output file path for saving plot
        log_scale: Use log scale for both axes
    """
    
    # Read CSV file
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: CSV file '{csv_file}' not found.")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: CSV file '{csv_file}' is empty.")
        return
    
    # Check if required columns exist
    if 'epsilon' not in df.columns:
        print("Error: 'epsilon' column not found in CSV file.")
        return
    
    # Find error column (could be L1_error, L2_error, or Linf_error)
    error_col = None
    for col in df.columns:
        if col.endswith('_error'):
            error_col = col
            break
    
    if error_col is None:
        print("Error: No error column found in CSV file.")
        return
    
    # Extract data
    epsilon = df['epsilon'].values
    error = df[error_col].values
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    if log_scale:
        plt.loglog(epsilon, error, 'bo-', linewidth=2, markersize=6)
        plt.xlabel('Epsilon (log scale)')
        plt.ylabel('Error (log scale)')
        plt.title('Gradient Calculation Error vs Epsilon Parameter (Log-Log Scale)')
    else:
        plt.plot(epsilon, error, 'bo-', linewidth=2, markersize=6)
        plt.xlabel('Epsilon')
        plt.ylabel('Error')
        plt.title('Gradient Calculation Error vs Epsilon Parameter')
    
    plt.grid(True, alpha=0.3)
    
    # Find optimal epsilon (minimum error)
    min_error_idx = np.argmin(error)
    optimal_epsilon = epsilon[min_error_idx]
    min_error = error[min_error_idx]
    
    # Highlight optimal point
    plt.plot(optimal_epsilon, min_error, 'ro', markersize=10, label=f'Optimal: ε={optimal_epsilon:.4f}, Error={min_error:.2f}')
    
    plt.legend()
    plt.tight_layout()
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    else:
        plt.show()
    
    # Print statistics
    print(f"\nEpsilon Error Analysis Statistics:")
    print(f"Data points: {len(epsilon)}")
    print(f"Epsilon range: [{epsilon[0]:.6f}, {epsilon[-1]:.6f}]")
    print(f"Error range: [{error.min():.2e}, {error.max():.2e}]")
    print(f"Optimal epsilon: {optimal_epsilon:.6f} (Error: {min_error:.2e})")

def main():
    parser = argparse.ArgumentParser(description='Plot epsilon error analysis results')
    parser.add_argument('csv_file', nargs='?', default='out/epsilon_error_analysis.csv',
                        help='Path to CSV file (default: out/epsilon_error_analysis.csv)')
    parser.add_argument('-o', '--output', help='Output file for plot (PNG format)')
    parser.add_argument('--linear', action='store_true', help='Use linear scale instead of log scale')
    
    args = parser.parse_args()
    
    # Check if CSV file exists
    if not os.path.exists(args.csv_file):
        print(f"Error: CSV file '{args.csv_file}' does not exist.")
        print("Please run the C++ program first to generate the data.")
        return 1
    
    # Plot the data
    plot_epsilon_error(args.csv_file, args.output, not args.linear)
    
    return 0

if __name__ == "__main__":
    exit(main())
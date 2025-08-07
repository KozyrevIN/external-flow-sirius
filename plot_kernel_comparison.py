#!/usr/bin/env python3
"""
Plot kernel comparison analysis with 2x2 subplot layout.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys

def plot_kernel_comparison(csv_files, output_file=None, function_name="f"):
    """
    Create 2x2 subplot comparing different kernels.
    
    Args:
        csv_files: List of CSV files for each kernel
        output_file: Output file path for saving plot
        function_name: Name of the function being analyzed
    """
    
    if len(csv_files) != 4:
        print(f"Error: Expected 4 CSV files for kernel comparison, got {len(csv_files)}")
        return
    
    kernel_names = ["kernel_2", "kernel_4", "kernel_6", "kernel_8"]
    
    # Read CSV files
    data = []
    for i, csv_file in enumerate(csv_files):
        try:
            df = pd.read_csv(csv_file)
            if 'epsilon' not in df.columns:
                print(f"Error: 'epsilon' column not found in {csv_file}")
                return
            
            # Find error column
            error_col = None
            for col in df.columns:
                if col.endswith('_error'):
                    error_col = col
                    break
            
            if error_col is None:
                print(f"Error: No error column found in {csv_file}")
                return
                
            data.append({
                'epsilon': df['epsilon'].values,
                'error': df[error_col].values,
                'kernel': kernel_names[i]
            })
            
        except FileNotFoundError:
            print(f"Error: CSV file '{csv_file}' not found.")
            return
        except pd.errors.EmptyDataError:
            print(f"Error: CSV file '{csv_file}' is empty.")
            return
    
    # Create 2x2 subplot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Kernel Comparison: {function_name} Function (L2 Norm)', fontsize=16, fontweight='bold')
    
    # Colors for each kernel
    colors = ['blue', 'red', 'green', 'orange']
    
    for i, (ax, kernel_data, color) in enumerate(zip(axes.flat, data, colors)):
        epsilon = kernel_data['epsilon']
        error = kernel_data['error']
        kernel_name = kernel_data['kernel']
        
        # Plot with log scale
        ax.loglog(epsilon, error, 'o-', color=color, linewidth=2, markersize=4)
        
        # Find optimal point
        min_error_idx = np.argmin(error)
        optimal_epsilon = epsilon[min_error_idx]
        min_error = error[min_error_idx]
        
        # Highlight optimal point
        ax.plot(optimal_epsilon, min_error, 'o', color='red', markersize=8,
                markeredgecolor='black', markeredgewidth=1)
        
        # Set labels and title
        ax.set_xlabel('Epsilon', fontsize=10)
        ax.set_ylabel('L2 Error', fontsize=10)
        ax.set_title(f'{kernel_name}\nOptimal: ε={optimal_epsilon:.4f}, Error={min_error:.2e}', 
                    fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Set consistent axis limits for comparison
        if i == 0:
            # Use first subplot to determine global limits
            global_epsilon_min = min([d['epsilon'].min() for d in data])
            global_epsilon_max = max([d['epsilon'].max() for d in data])
            global_error_min = min([d['error'].min() for d in data])
            global_error_max = max([d['error'].max() for d in data])
        
        ax.set_xlim(global_epsilon_min * 0.8, global_epsilon_max * 1.2)
        ax.set_ylim(global_error_min * 0.5, global_error_max * 2.0)
    
    plt.tight_layout()
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Kernel comparison plot saved to: {output_file}")
    else:
        plt.show()
    
    # Print comparison statistics
    print(f"\nKernel Comparison Statistics for {function_name}:")
    print("-" * 60)
    for i, kernel_data in enumerate(data):
        epsilon = kernel_data['epsilon']
        error = kernel_data['error']
        kernel_name = kernel_data['kernel']
        
        min_error_idx = np.argmin(error)
        optimal_epsilon = epsilon[min_error_idx]
        min_error = error[min_error_idx]
        
        print(f"{kernel_name:10} | Optimal ε: {optimal_epsilon:8.6f} | Min Error: {min_error:.3e}")

def main():
    parser = argparse.ArgumentParser(description='Plot kernel comparison analysis (2x2 layout)')
    parser.add_argument('csv_files', nargs=4, help='Four CSV files for kernel_2, kernel_4, kernel_6, kernel_8')
    parser.add_argument('-o', '--output', help='Output file for plot (PNG format)')
    parser.add_argument('--function', default='f', help='Function name for plot title')
    
    args = parser.parse_args()
    
    # Check if all CSV files exist
    for csv_file in args.csv_files:
        if not os.path.exists(csv_file):
            print(f"Error: CSV file '{csv_file}' does not exist.")
            return 1
    
    # Create the comparison plot
    plot_kernel_comparison(args.csv_files, args.output, args.function)
    
    return 0

if __name__ == "__main__":
    exit(main())
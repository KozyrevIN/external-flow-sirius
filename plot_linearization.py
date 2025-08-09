#!/usr/bin/env python3
"""
Plot epsilon error analysis with linear regression analysis.
Finds optimal epsilon and fits linear functions on both sides.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from scipy import stats

def linear_regression(x, y):
    """
    Perform linear regression on log-scale data.
    Returns slope, intercept, and r-squared value.
    """
    # Convert to log scale for linear regression
    log_x = np.log10(x)
    log_y = np.log10(y)
    
    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)
    r_squared = r_value ** 2
    
    return slope, intercept, r_squared

def plot_linearization_analysis(csv_file, output_file=None, eps_linear_min=None, eps_linear_max=None,
                               function_name="f", kernel_name="kernel", norm_name="L2"):
    """
    Plot epsilon vs error with linear regression analysis.
    
    Args:
        csv_file: Path to CSV file with epsilon and error columns
        output_file: Optional output file path for saving plot
        eps_linear_min: Minimum epsilon for linear regression range
        eps_linear_max: Maximum epsilon for linear regression range
        function_name: Name of the function being analyzed
        kernel_name: Name of the kernel used
        norm_name: Name of the norm used
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
    
    # Find error column
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
    
    # Find optimal epsilon (minimum error)
    min_error_idx = np.argmin(error)
    optimal_epsilon = epsilon[min_error_idx]
    min_error = error[min_error_idx]
    
    print(f"Optimal epsilon found: {optimal_epsilon:.6f} with {norm_name} error: {min_error:.3e}")
    
    # Set default linear regression ranges if not provided
    if eps_linear_min is None:
        eps_linear_min = epsilon[0]
    if eps_linear_max is None:
        eps_linear_max = epsilon[-1]
    
    # Split data into two regions for linear regression
    # Region 1: from eps_linear_min to optimal_epsilon
    mask1 = (epsilon >= eps_linear_min) & (epsilon <= optimal_epsilon)
    eps1 = epsilon[mask1]
    err1 = error[mask1]
    
    # Region 2: from optimal_epsilon to eps_linear_max
    mask2 = (epsilon >= optimal_epsilon) & (epsilon <= eps_linear_max)
    eps2 = epsilon[mask2]
    err2 = error[mask2]
    
    # Perform linear regression on both regions
    slope1, intercept1, r2_1 = None, None, None
    slope2, intercept2, r2_2 = None, None, None
    
    if len(eps1) >= 2:
        slope1, intercept1, r2_1 = linear_regression(eps1, err1)
        print(f"Region 1 linear fit: slope = {slope1:.3f}, R² = {r2_1:.3f}")
    
    if len(eps2) >= 2:
        slope2, intercept2, r2_2 = linear_regression(eps2, err2)
        print(f"Region 2 linear fit: slope = {slope2:.3f}, R² = {r2_2:.3f}")
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Plot original data
    plt.loglog(epsilon, error, 'bo-', linewidth=2, markersize=6, label='Data', alpha=0.8)
    
    # Highlight optimal point
    plt.plot(optimal_epsilon, min_error, 'ro', markersize=12, 
             label=f'Optimal: ε={optimal_epsilon:.6f}', markeredgecolor='black', markeredgewidth=2)
    
    # Plot linear regression lines (thick lines on top of data points)
    colors = ['green', 'purple']
    
    if slope1 is not None:
        # Generate points for first linear fit
        eps_fit1 = np.logspace(np.log10(eps_linear_min), np.log10(optimal_epsilon), 100)
        err_fit1 = 10**(intercept1 + slope1 * np.log10(eps_fit1))
        plt.plot(eps_fit1, err_fit1, '-', color=colors[0], linewidth=3, alpha=1.0, zorder=10,
                label=f'Linear fit 1: slope={slope1:.2f}')
    
    if slope2 is not None:
        # Generate points for second linear fit
        eps_fit2 = np.logspace(np.log10(optimal_epsilon), np.log10(eps_linear_max), 100)
        err_fit2 = 10**(intercept2 + slope2 * np.log10(eps_fit2))
        plt.plot(eps_fit2, err_fit2, '-', color=colors[1], linewidth=3, alpha=1.0, zorder=10,
                label=f'Linear fit 2: slope={slope2:.2f}')
    
    # Mark linear regression boundaries
    plt.axvline(x=eps_linear_min, color='gray', linestyle=':', alpha=0.7, label='Linear range')
    plt.axvline(x=eps_linear_max, color='gray', linestyle=':', alpha=0.7)
    
    plt.xlabel('Epsilon (log scale)', fontsize=12)
    plt.ylabel(f'{norm_name} Error (log scale)', fontsize=12)
    plt.title(f'Linearization Analysis: {function_name} with {kernel_name}', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)
    
    # Create summary text box
    summary_lines = [f"Optimal ε: {optimal_epsilon:.6f} ({norm_name} Error: {min_error:.3e})"]
    
    if slope1 is not None:
        summary_lines.append(f"Left slope: {slope1:.3f} (R²: {r2_1:.3f})")
    else:
        summary_lines.append("Left slope: N/A (insufficient data)")
    
    if slope2 is not None:
        summary_lines.append(f"Right slope: {slope2:.3f} (R²: {r2_2:.3f})")
    else:
        summary_lines.append("Right slope: N/A (insufficient data)")
    
    summary_text = " | ".join(summary_lines)
    
    plt.figtext(0.5, 0.02, summary_text, ha='center', va='bottom', fontsize=11,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)  # Make room for summary text
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Linearization plot saved to: {output_file}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot epsilon error analysis with linearization')
    parser.add_argument('csv_file', help='Path to CSV file with epsilon and error data')
    parser.add_argument('-o', '--output', help='Output file for plot (PNG format)')
    parser.add_argument('--eps-linear-min', type=float, help='Minimum epsilon for linear regression')
    parser.add_argument('--eps-linear-max', type=float, help='Maximum epsilon for linear regression')
    parser.add_argument('--function', default='f', help='Function name for plot title')
    parser.add_argument('--kernel', default='kernel', help='Kernel name for plot title')
    parser.add_argument('--norm', default='L2', help='Norm type for plot labels')
    
    args = parser.parse_args()
    
    # Check if CSV file exists
    if not os.path.exists(args.csv_file):
        print(f"Error: CSV file '{args.csv_file}' does not exist.")
        return 1
    
    # Create the linearization plot
    plot_linearization_analysis(args.csv_file, args.output, args.eps_linear_min, args.eps_linear_max,
                               args.function, args.kernel, args.norm)
    
    return 0

if __name__ == "__main__":
    exit(main())
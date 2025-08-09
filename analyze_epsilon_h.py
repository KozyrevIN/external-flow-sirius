#!/usr/bin/env python3
"""
Analysis script for epsilon vs h relationship from sphere error analysis.

This script analyzes the optimal epsilon values as a function of mesh size h
and provides various visualizations and trend analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

def load_data(csv_path):
    """Load the epsilon vs h analysis data."""
    try:
        df = pd.read_csv(csv_path)
        print(f"Loaded data with {len(df)} mesh sizes")
        print(f"h_max range: {df['h_max'].min():.4f} - {df['h_max'].max():.4f}")
        print(f"Optimal epsilon range: {df['optimal_epsilon'].min():.4f} - {df['optimal_epsilon'].max():.4f}")
        print(f"ε/h ratio range: {df['epsilon_h_ratio'].min():.2f} - {df['epsilon_h_ratio'].max():.2f}")
        return df
    except FileNotFoundError:
        print(f"Error: Could not find {csv_path}")
        print("Please run the sphere_error_analysis_example first to generate the data.")
        return None

def analyze_power_law_trend(df):
    """Analyze the power law trend in epsilon vs h data."""
    log_h = np.log10(df['h_max'])
    log_eps = np.log10(df['optimal_epsilon'])
    slope, intercept, r_value, p_value, std_err = linregress(log_h, log_eps)
    
    # Calculate coefficient errors
    n = len(df)
    slope_error = std_err
    intercept_error = std_err * np.sqrt(np.sum(log_h**2) / n)
    
    print(f"\n=== Power Law Trend Analysis ===")
    print(f"ε ∝ h^{slope:.4f} ± {slope_error:.4f}")
    print(f"Prefactor: 10^({intercept:.4f} ± {intercept_error:.4f}) = {10**intercept:.4f}")
    print(f"R² = {r_value**2:.4f}")
    print(f"p-value = {p_value:.2e}")
    
    return slope, intercept, r_value, slope_error, intercept_error

def plot_error_analysis(df):
    """Analyze the optimal errors achieved."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Error vs h (convergence analysis)
    ax1.loglog(df['h_max'], df['optimal_error'], 'o', markersize=6, alpha=0.7, color='purple', label='Optimal L2 Error')
    ax1.set_xlabel('h_max (log scale)')
    ax1.set_ylabel('Optimal L2 Error (log scale)')
    ax1.set_title('Error Convergence Analysis\n(Log-Log Scale)')
    ax1.grid(True, alpha=0.3)
    
    # Add convergence rate trend line
    log_h = np.log10(df['h_max'])
    log_err = np.log10(df['optimal_error'])
    slope_err, intercept_err, r_err, _, std_err_slope = linregress(log_h, log_err)
    
    h_trend = np.logspace(np.log10(df['h_max'].min()), np.log10(df['h_max'].max()), 100)
    err_trend = 10**(slope_err * np.log10(h_trend) + intercept_err)
    ax1.plot(h_trend, err_trend, '--', alpha=0.8, color='red', linewidth=2,
             label=f'Trend: Error ∝ h^{slope_err:.3f}±{std_err_slope:.3f} (R²={r_err**2:.3f})')
    ax1.legend()
    
    # Plot 2: Triangle count vs h
    ax2.loglog(df['h_max'], df['triangle_count'], '^', markersize=6, alpha=0.7, color='orange', label='Triangle Count')
    ax2.set_xlabel('h_max (log scale)')
    ax2.set_ylabel('Triangle Count (log scale)')
    ax2.set_title('Mesh Refinement\n(Log-Log Scale)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    
    # Save plot
    output_path = 'plots/error_convergence_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved error analysis to: {output_path}")
    
    print(f"\n=== Convergence Analysis ===")
    print(f"Error convergence rate: Error ∝ h^{slope_err:.4f} ± {std_err_slope:.4f}")
    print(f"R² = {r_err**2:.4f}")
    
    return slope_err, r_err

def fit_power_law(df):
    """Fit different functional forms to the epsilon-h relationship."""
    
    def power_law(h, a, b):
        """ε = a * h^b"""
        return a * (h ** b)
    
    def linear_law(h, a, b):
        """ε = a * h + b"""
        return a * h + b
    
    def sqrt_law(h, a, b):
        """ε = a * sqrt(h) + b"""
        return a * np.sqrt(h) + b
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Original data
    h_data = df['h_max'].values
    eps_data = df['optimal_epsilon'].values
    
    # Fit power law
    try:
        popt_power, pcov_power = curve_fit(power_law, h_data, eps_data, p0=[0.3, 0.5])
        a_power, b_power = popt_power
        a_power_err, b_power_err = np.sqrt(np.diag(pcov_power))
        
        # Fit linear law
        popt_linear, pcov_linear = curve_fit(linear_law, h_data, eps_data)
        a_linear, b_linear = popt_linear
        a_linear_err, b_linear_err = np.sqrt(np.diag(pcov_linear))
        
        # Fit sqrt law
        popt_sqrt, pcov_sqrt = curve_fit(sqrt_law, h_data, eps_data)
        a_sqrt, b_sqrt = popt_sqrt
        a_sqrt_err, b_sqrt_err = np.sqrt(np.diag(pcov_sqrt))
        
        # Generate smooth curves for plotting
        h_smooth = np.linspace(h_data.min(), h_data.max(), 200)
        
        # Plot 1: Linear scale comparison
        ax1.scatter(df['h_max'], df['optimal_epsilon'], s=40, alpha=0.7, label='Data', color='black', zorder=5)
        ax1.plot(h_smooth, power_law(h_smooth, a_power, b_power), '-', linewidth=2, 
                 label=f'Power law: ε = ({a_power:.3f}±{a_power_err:.3f})·h^({b_power:.3f}±{b_power_err:.3f})', alpha=0.8)
        ax1.plot(h_smooth, linear_law(h_smooth, a_linear, b_linear), '--', linewidth=2,
                 label=f'Linear: ε = ({a_linear:.3f}±{a_linear_err:.3f})·h + ({b_linear:.3f}±{b_linear_err:.3f})', alpha=0.8)
        ax1.plot(h_smooth, sqrt_law(h_smooth, a_sqrt, b_sqrt), '-.', linewidth=2,
                 label=f'Square root: ε = ({a_sqrt:.3f}±{a_sqrt_err:.3f})·√h + ({b_sqrt:.3f}±{b_sqrt_err:.3f})', alpha=0.8)
        ax1.set_xlabel('h_max')
        ax1.set_ylabel('Optimal ε')
        ax1.set_title('Functional Form Fitting\n(Linear Scale)')
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Log scale comparison
        ax2.loglog(df['h_max'], df['optimal_epsilon'], 'o', markersize=6, alpha=0.7, label='Data', color='black', zorder=5)
        ax2.loglog(h_smooth, power_law(h_smooth, a_power, b_power), '-', linewidth=2, 
                   label=f'Power law: ε = ({a_power:.3f}±{a_power_err:.3f})·h^({b_power:.3f}±{b_power_err:.3f})', alpha=0.8)
        ax2.set_xlabel('h_max (log scale)')
        ax2.set_ylabel('Optimal ε (log scale)')
        ax2.set_title('Power Law Fit\n(Log-Log Scale)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_path = 'plots/functional_form_fitting.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved functional form analysis to: {output_path}")
        
        # Calculate R² values
        eps_power_pred = power_law(h_data, a_power, b_power)
        eps_linear_pred = linear_law(h_data, a_linear, b_linear)
        eps_sqrt_pred = sqrt_law(h_data, a_sqrt, b_sqrt)
        
        ss_tot = np.sum((eps_data - np.mean(eps_data))**2)
        r2_power = 1 - np.sum((eps_data - eps_power_pred)**2) / ss_tot
        r2_linear = 1 - np.sum((eps_data - eps_linear_pred)**2) / ss_tot
        r2_sqrt = 1 - np.sum((eps_data - eps_sqrt_pred)**2) / ss_tot
        
        print(f"\n=== Functional Form Analysis ===")
        print(f"Power law:    ε = ({a_power:.4f} ± {a_power_err:.4f}) · h^({b_power:.4f} ± {b_power_err:.4f})    (R² = {r2_power:.4f})")
        print(f"Linear:       ε = ({a_linear:.4f} ± {a_linear_err:.4f}) · h + ({b_linear:.4f} ± {b_linear_err:.4f})   (R² = {r2_linear:.4f})")
        print(f"Square root:  ε = ({a_sqrt:.4f} ± {a_sqrt_err:.4f}) · √h + ({b_sqrt:.4f} ± {b_sqrt_err:.4f})    (R² = {r2_sqrt:.4f})")
        
        return {'power': (a_power, b_power, r2_power, a_power_err, b_power_err), 
                'linear': (a_linear, b_linear, r2_linear, a_linear_err, b_linear_err),
                'sqrt': (a_sqrt, b_sqrt, r2_sqrt, a_sqrt_err, b_sqrt_err)}
        
    except Exception as e:
        print(f"Error in curve fitting: {e}")
        return None

def main():
    """Main analysis function."""
    print("=== Epsilon vs H Analysis ===")
    
    # Ensure plots directory exists
    import os
    os.makedirs('plots', exist_ok=True)
    
    # Load data
    df = load_data('out/epsilon_vs_h_analysis.csv')
    if df is None:
        return
    
    print(f"\nDataset: {len(df)} mesh sizes, h ∈ [{df['h_max'].min():.4f}, {df['h_max'].max():.4f}]")
    print(f"Epsilon range: {df['optimal_epsilon'].min():.4f} - {df['optimal_epsilon'].max():.4f}")
    
    # 1. Power law trend analysis
    print(f"\n1. Analyzing power law trend...")
    slope, intercept, r_value, slope_err, intercept_err = analyze_power_law_trend(df)
    
    # 2. Error convergence analysis  
    print(f"\n2. Analyzing error convergence...")
    slope_err, r_err = plot_error_analysis(df)
    
    # 3. Functional form fitting
    print(f"\n3. Fitting functional forms...")
    fit_results = fit_power_law(df)
    
    # Summary
    print(f"\n=== SUMMARY ===")
    print(f"Dataset: {len(df)} mesh sizes, h ∈ [{df['h_max'].min():.4f}, {df['h_max'].max():.4f}]")
    print(f"Optimal epsilon relationship: ε ∝ h^{slope:.4f} ± {slope_err:.4f} (R² = {r_value**2:.4f})")
    print(f"Error convergence: Error ∝ h^{slope_err:.3f} (R² = {r_err**2:.4f})")
    
    if fit_results:
        best_fit = max(fit_results.items(), key=lambda x: x[1][2])
        fit_name, fit_params = best_fit
        if len(fit_params) >= 5:  # With error estimates
            a, b, r2, a_err, b_err = fit_params
            print(f"Best functional form: {fit_name} with R² = {r2:.4f}")
        else:
            a, b, r2 = fit_params[:3]
            print(f"Best functional form: {fit_name} with R² = {r2:.4f}")
    
    print(f"\nPlots saved to plots/ directory")

if __name__ == "__main__":
    main()
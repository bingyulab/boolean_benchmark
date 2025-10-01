import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import argparse
import numpy as np

# Set up the plotting style
plt.style.use('seaborn-v0_8')

# Define constants
METHODS = ['CASPO', 'VNS', 'GA', 'ILP']
COLORS = ['#E41A1C', '#377EB8', '#4DAF4A', '#FF7F00']
METHOD_COLORS = dict(zip(METHODS, COLORS))

METRICS_CONFIG = {
    'attractors': {
        'metrics': ['jaccard', 'hamming'],
        'titles': ['Jaccard Similarity', 'Hamming Similarity'],
        'markers': ['o', 's']
    },
    'coverage': {
        'metrics': ['f1_score', 'mse'],
        'titles': ['F1 Score', 'MSE'],
        'markers': ['^', 'D']
    },
    'topology': {
        'metrics': ['node_jaccard_topology', 'edge_jaccard_topology'],
        'titles': ['Node Jaccard', 'Edge Jaccard'],
        'markers': ['o', 's']
    },
    'performance': {
        'metrics': ['total_time', 'recon_total'],
        'titles': ['Runtime (log scale)', 'Reconstructed Attractors'],
        'markers': ['o', 's']
    }
}

AVG_COLUMNS = [
    'jaccard', 'hamming', 'f1_score', 'mse',
    'recon_total', 'total_time', 
    'node_jaccard_topology', 'edge_jaccard_topology'
]

GROUP_COLUMNS = ['method', 'change_percent']

def load_and_aggregate_data(pattern: str, param_label: str, time_threshold: int = 500) -> pd.DataFrame:
    """Load and aggregate data for a single parameter set."""
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No files match pattern {pattern}")
    
    # Load and concatenate all files
    df = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    
    # Filter by time threshold and add parameter label
    df = df[df['total_time'] < time_threshold].copy()
    df['parameter_type'] = param_label
    
    # Group by method and change_percent, then calculate means
    return df.groupby(GROUP_COLUMNS + ['parameter_type'], dropna=False)[AVG_COLUMNS].mean().reset_index()

def load_and_process_data(pattern_50: str, pattern_240: str, time_threshold: int = 500) -> pd.DataFrame:
    """Load and combine data from two parameter sets with optimized processing."""
    df_50 = load_and_aggregate_data(pattern_50, 'parameter_50', time_threshold)
    df_240 = load_and_aggregate_data(pattern_240, 'parameter_240', time_threshold)
    
    return pd.concat([df_50, df_240], ignore_index=True)

def plot_metric_comparison(df: pd.DataFrame, ax, metric: str, title: str, marker: str = 'o'):
    """Optimized plotting function with vectorized operations."""
    # Pre-filter data for efficiency
    metric_data = df[['method', 'parameter_type', 'change_percent', metric]].dropna()
    
    for method in METHODS:
        method_data = metric_data[metric_data['method'] == method]
        if method_data.empty:
            continue
            
        for param_type, linestyle, param_label in [
            ('parameter_50', '-', '50'),
            ('parameter_240', '--', '240')
        ]:
            param_data = method_data[method_data['parameter_type'] == param_type]
            if not param_data.empty:
                # Sort by change_percent for smooth lines
                param_data = param_data.sort_values('change_percent')
                ax.plot(
                    param_data['change_percent'], param_data[metric],
                    marker=marker, linewidth=1.5, markersize=3,
                    label=f"{method} ({param_label})",
                    color=METHOD_COLORS[method], linestyle=linestyle
                )
    
    ax.set_xlabel('Change Percentage')
    ax.set_ylabel(metric.replace('_', ' ').title())
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

def create_comparison_plots(df: pd.DataFrame, dataset_name: str, output_dir: Path):
    """Create all comparison plots using configuration-driven approach."""
    
    for plot_type, config in METRICS_CONFIG.items():
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f'{plot_type.title()} Metrics - Parameter Comparison ({dataset_name})', 
                    fontsize=16, fontweight='bold')
        
        for i, (metric, title, marker) in enumerate(zip(config['metrics'], config['titles'], config['markers'])):
            plot_metric_comparison(df, axes[i], metric, title, marker)
            
            # Special handling for log scale
            if metric == 'total_time':
                axes[i].set_yscale('log')
        
        plt.tight_layout()
        fig.savefig(output_dir / f'parameter_comparison_{plot_type}.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

def create_summary_comparison(df: pd.DataFrame, dataset_name: str, output_dir: Path) -> pd.DataFrame:
    """Create optimized summary statistics comparison."""
    # Use groupby for efficient aggregation
    summary_df = df.groupby(['parameter_type', 'method'])[
        ['jaccard', 'f1_score', 'total_time', 'recon_total']
    ].mean().reset_index()
    
    # Rename columns for clarity
    summary_df.rename(columns={
        'jaccard': 'avg_jaccard',
        'f1_score': 'avg_f1_score', 
        'total_time': 'avg_runtime',
        'recon_total': 'avg_recon_total'
    }, inplace=True)
    
    summary_df.to_csv(output_dir / 'parameter_comparison_summary.csv', index=False)
    return summary_df

def create_statistical_comparison(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """Create detailed statistical comparison between parameter sets."""
    stats_comparison = []
    
    for method in METHODS:
        method_data = df[df['method'] == method]
        if method_data.empty:
            continue
            
        param_50 = method_data[method_data['parameter_type'] == 'parameter_50']
        param_240 = method_data[method_data['parameter_type'] == 'parameter_240']
        
        for metric in ['jaccard', 'f1_score', 'total_time', 'recon_total']:
            if not param_50.empty and not param_240.empty:
                stats = {
                    'method': method,
                    'metric': metric,
                    'param_50_mean': param_50[metric].mean(),
                    'param_240_mean': param_240[metric].mean(),
                    'param_50_std': param_50[metric].std(),
                    'param_240_std': param_240[metric].std(),
                    'improvement_ratio': param_240[metric].mean() / param_50[metric].mean() if param_50[metric].mean() != 0 else np.nan
                }
                stats_comparison.append(stats)
    
    stats_df = pd.DataFrame(stats_comparison)
    stats_df.to_csv(output_dir / 'parameter_statistical_comparison.csv', index=False)
    return stats_df

def main():
    parser = argparse.ArgumentParser(description='Compare results between two parameter sets')
    parser.add_argument('--dataset', default='TCell', help='Dataset name')
    parser.add_argument('--pattern_50', default='output/comparison_TCell_*.csv', 
                       help='Pattern for parameter 50 files')
    parser.add_argument('--pattern_240', default='output/candidate/comparison_TCell_*.csv',
                       help='Pattern for parameter 240 files')
    parser.add_argument('--output_dir', default='output/parameters', 
                       help='Output directory')
    parser.add_argument('--time_threshold', type=int, default=500,
                       help='Time threshold for filtering results')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Load and process data
        print("Loading and processing data...")
        df = load_and_process_data(args.pattern_50, args.pattern_240, args.time_threshold)
        
        # Create plots
        print("Creating comparison plots...")
        create_comparison_plots(df, args.dataset, output_dir)
        
        # Create summaries
        print("Creating summary statistics...")
        summary_df = create_summary_comparison(df, args.dataset, output_dir)
        stats_df = create_statistical_comparison(df, output_dir)
        
        print(f"\nParameter comparison analysis completed!")
        print(f"Results saved to '{output_dir}'")
        print("\nGenerated files:")
        for plot_type in METRICS_CONFIG.keys():
            print(f"  - parameter_comparison_{plot_type}.png")
        print("  - parameter_comparison_summary.csv")
        print("  - parameter_statistical_comparison.csv")
        
        # Print quick summary
        print(f"\nDataset: {args.dataset}")
        print(f"Total records processed: {len(df)}")
        print(f"Methods analyzed: {df['method'].nunique()}")
        print(f"Parameter types: {', '.join(df['parameter_type'].unique())}")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
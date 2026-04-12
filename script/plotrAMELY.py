import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

def process_files(stats_file, female_max_rAMELY, male_min_rAMELY, lowest_unique_peptides_num):

    df = pd.read_csv(stats_file, dtype={
        'Individual': str,
        'Mean': float,
        'unqiue_AMELYX_peptides': float,
        'total_amelogenin_peptides': float
    })
    
    df['Sex'] = 'U'
    df.loc[df['Mean'] <= female_max_rAMELY, 'Sex'] = 'F'
    df.loc[df['Mean'] >= male_min_rAMELY, 'Sex'] = 'M'
    df.loc[(df['Mean'] <= female_max_rAMELY) & (df['unqiue_AMELYX_peptides'] < lowest_unique_peptides_num), 'Sex'] = 'F(High risk)'
    df.loc[(df['Mean'] >= male_min_rAMELY) & (df['unqiue_AMELYX_peptides'] < lowest_unique_peptides_num), 'Sex'] = 'M(High risk)'

    df.columns = ['Individual', 'rAMELY_mean', 'unqiue_AMELYX_peptides', 'total_amelogenin_peptides', 'Sex']
    new_column_order = ['Individual', 'Sex', 'rAMELY_mean', 'unqiue_AMELYX_peptides', 'total_amelogenin_peptides']
    df_reordered = df[new_column_order]
    df_reordered.to_csv("Sex_assessment_report.csv", index=False)

    return None

def calculate_figsize(num_individuals):
    base_height_per_row = 0.4  
    height = max(6, min(20, num_individuals * base_height_per_row + 2))
    width = 6  
    return (width, height)

def rAMELYPlot(csv_file, female_max_rAMELY, male_min_rAMELY, output_plot=None, dpi=300, use_seaborn_style=True):

    # Read data
    df = pd.read_csv(csv_file)

    # Check required columns
    required_columns = ['Individual', 'Sex', 'rAMELY_mean', 'unqiue_AMELYX_peptides', 'total_amelogenin_peptides']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    figsize = calculate_figsize(len(df))

    # Set style if requested
    if use_seaborn_style:
        sns.set_style("whitegrid")
    
    # Create figure
    _fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Get unique sexes for coloring
    face_color_map = {'F': '#377EB8', 'M': '#E41A1C', 'F(High risk)': '#7F7F7F','M(High risk)': '#7F7F7F','U': '#7F7F7F'}
    edge_color_map = {'F': '#377EB8', 'M': '#E41A1C', 'F(High risk)': '#377EB8','M(High risk)': '#E41A1C','U': '#7F7F7F'}

    # Plot each point with error bars
    for i, (idx, row) in enumerate(df.iterrows()):
        ax.errorbar(
            x=row['rAMELY_mean'],
            y=i,
            fmt='o',
            markerfacecolor=face_color_map[row['Sex']],  # Face color
            markeredgecolor=edge_color_map[row['Sex']],  # Edge color
            markersize=8,
            markeredgewidth=2,
            capsize=5,
            label=row['Sex'] if row['Sex'] not in df.iloc[:i]['Sex'].values else ""
        )
    
    # Show unique peptides number in right axis
    secax = ax.secondary_yaxis('right')
    y_labels = [f"n={int(row['unqiue_AMELYX_peptides'])}" for _, row in df.iterrows()]
    secax.set_yticks(range(len(df)))
    secax.set_yticklabels(y_labels)

    # Add vertical dashed lines
    ax.axvline(x=female_max_rAMELY, color=face_color_map['F'], linestyle='--', alpha=0.7, linewidth=1.5)
    ax.axvline(x=male_min_rAMELY, color=face_color_map['M'], linestyle='--', alpha=0.7, linewidth=1.5)

    ax.set_ylim(-0.5, len(df) - 0.5)  # Remove extra whitespace

    # Get limit of axes
    x_min, x_max = ax.get_xlim()
    _y_min, y_max = ax.get_ylim()

    # Add tag
    ax.text((x_min + female_max_rAMELY)/2, y_max * 1.001, 'XX(Female)', 
        ha='center', va='bottom', fontsize=12)
    ax.text((male_min_rAMELY + x_max)/2, y_max * 1.001, 'XY(Male)', 
        ha='center', va='bottom', fontsize=12)
    
    # Add legend
    legend_elements = [
    Line2D([0], [0], marker='o', color='w', 
           markerfacecolor=face_color_map['F'], 
           markeredgecolor=edge_color_map['F'],
           markersize=8, markeredgewidth=2, label='F'),
    Line2D([0], [0], marker='o', color='w', 
           markerfacecolor=face_color_map['F(High risk)'], 
           markeredgecolor=edge_color_map['F(High risk)'],
           markersize=8, markeredgewidth=2, label='F(High risk)'),
    Line2D([0], [0], marker='o', color='w', 
           markerfacecolor=face_color_map['M'], 
           markeredgecolor=edge_color_map['M'],
           markersize=8, markeredgewidth=2, label='M'),
    Line2D([0], [0], marker='o', color='w', 
           markerfacecolor=face_color_map['M(High risk)'], 
           markeredgecolor=edge_color_map['M(High risk)'],
           markersize=8, markeredgewidth=2, label='M(High risk)'),
    Line2D([0], [0], marker='o', color='w', 
           markerfacecolor=face_color_map['U'], 
           markeredgecolor=edge_color_map['U'],
           markersize=8, markeredgewidth=2, label='U'),
    ]
    ax.legend(handles=legend_elements, title='Sex', loc='best')

    # Customize plot
    ax.set_xlabel(r'$R_{AMELY}$', fontsize=12, fontweight='bold')
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['Individual'])
    
    ax.grid(None)
    
    plt.tight_layout()
    
    if output_plot:
        plt.savefig(output_plot, bbox_inches='tight', dpi=dpi)
        print(f"rAMELY distribution plot saved to: {output_plot}")
    else:
        plt.show()
    
    plt.close()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='Process gender and rAMELY statistics files to calculate minimum lower bound of 95% CI for male individuals'
    )
    parser.add_argument('stats_file', type=str, help='Path to CSV file with rAMELY statistics')
    parser.add_argument('female_max_rAMELY', type=float, help='The highest rAMELY value for female')
    parser.add_argument('male_min_rAMELY', type=float, help='The lowest rAMELY value for male')
    parser.add_argument('minimum_unique_peptides_num', type=float, help='The minimum number of unique peptides required for a reliable rAMELY estimate')    
    args = parser.parse_args()

    stats_file = args.stats_file
    female_max_rAMELY = args.female_max_rAMELY
    male_min_rAMELY = args.male_min_rAMELY
    minimum_unique_peptides_num = args.minimum_unique_peptides_num
    
    process_files(stats_file, female_max_rAMELY, male_min_rAMELY, minimum_unique_peptides_num)

    rAMELYPlot("Sex_assessment_report.csv", female_max_rAMELY, male_min_rAMELY, 'rAMELY_SexDeterminer.pdf')
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def process_files(gender_file, stats_file):

    gender_df = pd.read_csv(gender_file, header=None)
    
    individuals = gender_df.iloc[:, 0].tolist()  
    genders = gender_df.iloc[:, 2].tolist()      
    
    gender_dict = {}
    for ind, gender in zip(individuals, genders):
        gender_dict[str(ind).strip()] = str(gender).strip()
    
    stats_df = pd.read_csv(stats_file)
    
    stats_df['Sex'] = stats_df.iloc[:, 0].apply(
        lambda x: gender_dict.get(str(x).strip())
    )
    
    male_stats = stats_df[stats_df['Sex'] == 'M']
    female_stats = stats_df[stats_df['Sex'] == 'F']

    male_min_rAMELY = male_stats['Mean'].min()
    female_max_rAMELY = female_stats['Mean'].max()
    lowest_unique_peptides_num = stats_df['unqiue_AMELYX_peptides'].min()

    stats_df.to_csv(stats_file, index=False)
    
    with open('rAMELY_threshold.txt', 'w') as f:
        f.write(f"The highest rAMLEY value for female is {female_max_rAMELY}\n")
        f.write(f"The lowest rAMLEY value for male is {male_min_rAMELY}\n")
        if female_max_rAMELY >= male_min_rAMELY:  
            f.write("!! WORNING: The maximum rAMELY value for female is greater than or equal to the minimum rAMELY value!")
            f.write("!! Please check the input data and re-evaluate the threshold!\n")
        f.write(f"The lowest unique AMELYX peptide number is {lowest_unique_peptides_num}\n")
    print(f"rAMELY threshold saved to: rAMELY_threshold.txt")

    return female_max_rAMELY, male_min_rAMELY

def calculate_figsize(num_individuals):
    base_height_per_row = 0.4  
    height = max(6, min(20, num_individuals * base_height_per_row + 2))
    width = 6  
    return (width, height)

def rAMELYPlot(csv_file, female_max_rAMELY, male_min_rAMELY, output_plot=None, dpi=300, use_seaborn_style=True):

    # Read data
    df = pd.read_csv(csv_file)

    # Check required columns
    required_columns = ['Individual', 'Mean', 'Sex', 'unqiue_AMELYX_peptides', 'total_amelogenin_peptides']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")
    
    figsize = calculate_figsize(len(df))

    # Set style if requested
    if use_seaborn_style:
        sns.set_style("whitegrid")
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Get unique sexes for coloring
    color_map = {'F': '#377EB8', 'M': '#E41A1C'}

    # Plot each point with error barsx
    for i, (idx, row) in enumerate(df.iterrows()):
        ax.errorbar(
            x=row['Mean'],
            y=i,
            fmt='o',
            color=color_map[row['Sex']],
            markersize=8,
            capsize=5,
            label=row['Sex'] if row['Sex'] not in df.iloc[:i]['Sex'].values else ""
        )
    
    # Show unique peptides number in right axis
    secax = ax.secondary_yaxis('right')
    y_labels = [f"n={int(row['unqiue_AMELYX_peptides'])}" for _, row in df.iterrows()]
    secax.set_yticks(range(len(df)))
    secax.set_yticklabels(y_labels)   

    # Add vertical dashed lines
    ax.axvline(x=female_max_rAMELY, color=color_map['F'], linestyle='--', alpha=0.7, linewidth=1.5)
    ax.axvline(x=male_min_rAMELY, color=color_map['M'], linestyle='--', alpha=0.7, linewidth=1.5)

    ax.set_ylim(-0.5, len(df) - 0.5)  # Remove extra whitespace

    # Get limit of axes
    x_min, x_max = ax.get_xlim()
    _y_min, y_max = ax.get_ylim()

    # Add tag
    ax.text((x_min + female_max_rAMELY)/2, y_max * 1.001, 'XX(Female)', 
        ha='center', va='bottom', fontsize=12)
    ax.text((male_min_rAMELY + x_max)/2, y_max * 1.001, 'XY(Male)', 
        ha='center', va='bottom', fontsize=12)

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
        description='Process gender and rAMELY statistics files to calculate the threshold for distinguishing male and female individuals'
    )
    parser.add_argument('gender_file', type=str, help='Path to CSV file with gender information (first row: individual names, second row: path of raw data, third row: genders)')
    parser.add_argument('stats_file', type=str, help='Path to CSV file with rAMELY statistics')
    args = parser.parse_args()

    gender_file = args.gender_file
    stats_file = args.stats_file
    
    female_max_rAMELY, male_min_rAMELY = process_files(gender_file, stats_file)

    rAMELYPlot(stats_file, female_max_rAMELY, male_min_rAMELY, 'rAMELY_Distribution.pdf')
#!/usr/bin/env python3 

import pandas as pd
import argparse

def rAMELY_computer(ind_name, df_fn):
    df = pd.read_csv(df_fn)
    print(f"Calculate the rAMELY ratio for {ind_name}")
    summary_df = df.groupby('Individual').agg(
        AMELY_unique_peptides=('Classification', lambda x: (x == 'AMELY-unique').sum()),
        AMELX_unique_peptides=('Classification', lambda x: (x == 'AMELX-unique').sum()),
        Both=('Classification', lambda x: (x == 'Both').sum()),
    )
    unique_peptides = (
        summary_df['AMELX_unique_peptides'].item() + 
        summary_df['AMELY_unique_peptides'].item()
    )
    total_peptides = (
        summary_df['AMELX_unique_peptides'].item() + 
        summary_df['AMELY_unique_peptides'].item() +
        summary_df['Both'].item()        
    )
    if summary_df['AMELY_unique_peptides'].item() > 0:
        summary_df['Mean'] = summary_df['AMELY_unique_peptides'] / unique_peptides
    else:
        summary_df['Mean'] = 0
    summary_df['unqiue_AMELYX_peptides'] = unique_peptides
    summary_df['total_amelogenin_peptides'] = total_peptides
    
    summary_df = summary_df.reset_index()
    res_df = summary_df[['Individual', 'Mean', 'unqiue_AMELYX_peptides', 'total_amelogenin_peptides']]
    res_df.to_csv(f'{ind_name}-rAMELY_summary_statistics.csv', index=False)        
    print(f"Results saved to: {ind_name}-rAMELY_summary_statistics.csv")
    
if __name__ == "__main__":    
    parser = argparse.ArgumentParser(
        description="Calculate the rAMELY ratio",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("ind_name", help="The name of sample")
    parser.add_argument("classified_peptides", help="The path of classified_peptides.csv")
    args = parser.parse_args()
    
    rAMELY_computer(args.ind_name, args.classified_peptides)
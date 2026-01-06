#!/usr/bin/env python
"""
Classify the peptides to `AMELY-unique`, `AMELX-unique`, or `Both` (means it matchs both AMELY and AMELX) categories
"""

import pandas as pd
import re
import numpy as np
import argparse
from Bio import SeqIO

def parse_ref_sequences(ref_fn: str) -> tuple:
    """Parse reference fasta file to get AMELX and AMELY peptide sequences"""
    AMELX_seqlist = []
    AMELY_seqlist = []
    with open(ref_fn, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            header = record.id + " " + record.description
            sequence = str(record.seq)
            gene = extract_gene_from_accession(header)
            if gene == 'AMELX':
                AMELX_seqlist.append(sequence)
            elif gene == 'AMELY':
                AMELY_seqlist.append(sequence)
    return AMELX_seqlist, AMELY_seqlist

def is_in_AMELXY_seqlist(peptide: str, AMELXY_seqlist: list) -> bool:
    """Check if peptide is in AMELX or AMELY sequence list"""
    for seq in AMELXY_seqlist:
        if seq.find(peptide) != -1:
            return True
    return False

def preprocess_peptide(peptide: str) -> str:
    """
    Process peptide sequence:
    1. If flanked by dots at second positions (e.g., .ABC.), extract middle portion
    2. Remove all parentheses and their contents (e.g., A(12.34)B → AB)
    """
    if len(peptide) >= 4 and peptide[1] == '.' and peptide[-2] == '.':
        peptide = peptide[2:-2]
    elif peptide[1] == '.' and not peptide[-2] == '.':
        peptide = peptide[2:]
    elif not peptide[1] == '.' and peptide[-2] == '.':
        peptide = peptide[:-2]
    peptide = re.sub(r'\([^)]*\)', '', peptide)
    return peptide

def extract_gene_from_accession(accession: str) -> str:
    """Extract gene name from Accession field when Gene column is empty"""
    # First, search for AMELX and AMELY directly
    if accession.find('AMELX') != -1:
        return 'AMELX'
    if accession.find('AMELY') != -1:
        return 'AMELY'

    # Try GN=([\w-]+) pattern 
    gn_match = re.search(r'.+\sGN=([\w-]+)\s*', accession)
    if gn_match:
        return gn_match.group(1)
    
    # Try ^>{0,1}([a-zA-Z0-9]+)(_|-)[\w-]+ prefix pattern
    prefix_match = re.search(r'^>{0,1}([a-zA-Z0-9]+)(_|-)[\w-]+', accession)
    if prefix_match:
        return prefix_match.group(1)

    return 'Unknown'

def process_file(file_path: str) -> pd.DataFrame:
    """
    Process a single search results file with cleaning and preprocessing steps.
    """
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()
    
    # Clean and preprocess data
    # df = df.drop_duplicates(subset=['Protein Group', '-10LgP'])
    df = df[~df['Accession'].astype(str).str.startswith('#CONTAM')]
    df = df[~df['Peptide'].astype(str).str.contains(r'\(ins\)|\(del\)', regex=True)]
    df['Peptide'] = df['Peptide'].apply(preprocess_peptide)
    if 'Gene' not in df.columns:
        df['Gene'] = np.nan
    df['Gene'] = np.where(
        df['Gene'].notna(),
        df['Gene'].astype(str).str.strip(),
        df['Accession'].apply(extract_gene_from_accession)
    )
    
    return df

def classify_peptides(df: pd.DataFrame, AMELX_seqlist: list, AMELY_seqlist: list) -> pd.DataFrame:
    """
    Classify peptides into AMELX-unique, AMELY-unique, or both (means it matchs both AMELY and AMELX) categories.
    Returns filtered DataFrame and weighted statistics.
    """
    # Get Spec and Intensity columns (adjust if multiple exist)
    if '#Spec' in df.columns:
        spec_cols = ['#Spec']
    else:
        spec_cols = [col for col in df.columns if col.startswith('#Spec ')]
    if 'Intensity' in df.columns:
        intensity_cols = ['Intensity']
    else:
        intensity_cols = [col for col in df.columns if col.startswith('Intensity ')]
    if 'Area' in df.columns:
        area_cols = ['Area']
    else:
        area_cols = [col for col in df.columns if col.startswith('Area ')]

    classified_rows = []
    
    for peptide, group in df.groupby('Peptide'):
        genes = set(group['Gene'].unique())
        classification = None
        
        if genes == {'AMELX'}:
            classification = 'AMELX-unique'
        elif genes == {'AMELY'}:
            classification = 'AMELY-unique'
        elif genes == {'AMELX', 'AMELY'}:
            classification = 'Both'
        else:
            continue
        
        if classification == 'AMELX-unique' and is_in_AMELXY_seqlist(peptide, AMELY_seqlist):
            classification = 'AMELX-unique-match-AMELY'
        if classification == 'AMELY-unique' and is_in_AMELXY_seqlist(peptide, AMELX_seqlist):
            classification = 'AMELY-unique-match-AMELX'

        # Take the first row for each peptide (already deduplicated)
        row = group.iloc[0].copy()
        columns_to_extract = ['Peptide', 'RT', 'Source File', 'Scan']
        existing_columns = [col for col in columns_to_extract if col in group.columns]
        row = row[existing_columns]
        row['Classification'] = classification
        
        if spec_cols:
            spec_sum = 0
            for spec_col in spec_cols:
                unique_values = group[spec_col].dropna().unique()
                spec_sum += unique_values.sum()
            row['#Spec'] = spec_sum
        if intensity_cols:
            intensity_sum = 0
            for intensity_col in intensity_cols:
                unique_values = group[intensity_col].replace({0: 1, np.nan: 1}).dropna().unique()
                intensity_sum += unique_values.sum()
            row['Intensity'] = intensity_sum
        if area_cols:
            area_sum = 0
            for area_col in area_cols:
                unique_values = group[area_col].replace({0: 1, np.nan: 1}).dropna().unique()
                area_sum += unique_values.sum()
            row['Area'] = area_sum
        # Add to classified results
        classified_rows.append(row)

    if not classified_rows:
        return pd.DataFrame({'Peptide': [], 'RT': [], 'Source File': [], 'Scan': [], 'Classification': []})
    
    classified_df = pd.DataFrame(classified_rows)
    return classified_df

def process_individual(ind_name: str, peptide_results: str, ref_fn: str) -> None:

    print(f"Classify the peptides of {ind_name} into AMELY-unique, AMELX-unqiue, and 'Both' category")    
    # Process and classify peptides
    df = process_file(peptide_results)
    
    AMELX_seqlist = []
    AMELY_seqlist = []

    AMELX_seqlist, AMELY_seqlist = parse_ref_sequences(ref_fn)

    classified_df = classify_peptides(df, AMELX_seqlist, AMELY_seqlist)
        
    # Add individual identifier and save classified peptides
    classified_df.insert(0, 'Individual', ind_name)
    
    # Save results if data exists
    classified_df.to_csv(f"{ind_name}-classified_peptides.csv", index=False)
    print(f"Classified peptides result saved to: {ind_name}-classified_peptides.csv")


if __name__ == "__main__":    
    parser = argparse.ArgumentParser(
        description="Process proteomics database searching data and classify the peptides",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("ind_name", help="The name of sample")
    parser.add_argument("peptides_result", help="The path of peptides database searching result file")
    parser.add_argument("ref_fn", help="The path of reference fasta file")
    args = parser.parse_args()
    
    process_individual(args.ind_name, args.peptides_result, args.ref_fn)
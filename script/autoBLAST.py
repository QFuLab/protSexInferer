#!/usr/bin/env python3

import pandas as pd
import argparse
from typing import Tuple, List
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import subprocess
import sys

def is_perfect_match(hsp, query_length: int) -> bool:
    if hsp.align_length != query_length:
        return False
        
    identity_percentage = (hsp.identities / hsp.align_length) * 100
    if identity_percentage < 100.0:
        return False
        
    if hsp.gaps > 0:
        return False
        
    return True

def BLAST_peptides(ind_name: str, temp_fasta_filename: str) -> List[Tuple[str, List]]:

    print(f"Run BLAST for {ind_name}...")
    output_xml = f"{ind_name}-BlastOutPut.xml"
    blastp_cmd = [
        "blastp",
        "-task", "blastp-short",           
        "-query", temp_fasta_filename,
        "-db", "nr",
        "-out", output_xml,
        "-outfmt", "5",     
    ]    
    
    print("BLAST command: ", " ".join(blastp_cmd))

    try:
        result = subprocess.run(
            blastp_cmd,
            check=True,      
            capture_output=True,  
            text=True        
    )
        print(f"Finish BLASTP of {ind_name}")
    
    except subprocess.CalledProcessError as e:
        print(f"Error！Exit Code: {e.returncode}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error：Can not found BLASTP program. Please install it on the PATH directory")
        sys.exit(1)

    results = []
    print(f"Parse BLAST results of {ind_name}...")
    with open(output_xml, 'r') as handle:
        blast_records = NCBIXML.parse(handle)
        for blast_record in blast_records:
            query = blast_record.query
            query_length = len(query)
            hits = []
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if is_perfect_match(hsp, query_length):
                        result = {
                            'hit_id': alignment.hit_id,
                            'hit_def': alignment.hit_def,
                            'evalue': hsp.expect,
                            'alignment_length': hsp.align_length,
                            'query': hsp.query,
                            'match': hsp.match,
                            'subject': hsp.sbjct
                        }
                        hits.append(result)
            results.append((query, hits))
    os.remove(output_xml)
    print(f"Finish BLAST parsing for {ind_name}")

    return results

def autoBLAST(ind_name: str, peptides: List[str]) -> List[Tuple[str, List]]:

    fn = f"{ind_name}-peptides-need-BLAST.fasta"
    records = []
    for peptide in peptides:
        my_seq = Seq(peptide)
        record = SeqRecord(my_seq, id=peptide, description="")
        records.append(record) 
    SeqIO.write(records, fn, "fasta")    
    
    results = BLAST_peptides(ind_name, fn)

    return results

def not_Amelogenin(hit_def: str) -> bool:
    hit_def = hit_def.lower()
    return (hit_def.find('hypothetical') != -1 and hit_def.find('predicted') != -1 and 
            hit_def.find('amelogenin') != -1 and hit_def.find('amely') != -1 and hit_def.find('amelx') != -1)

def get_peptides_to_remove(ind_name: str, results: List[Tuple[str, List]]) -> List:
    peptides_to_remove = []
    done = False

    for (peptide, hits) in results:
        for hit in hits:
            if not_Amelogenin(hit['hit_def']):
                peptides_to_remove.append(peptide)
                done = True
                break
            if done:
                break
                
    with open(f"{ind_name}-removed_peptides.list", "w") as f:
        f.write('\n'.join(peptides_to_remove))

    return peptides_to_remove

def output_BLAST_report(ind_name: str, results: List[Tuple[str, List]]):
    with open(f'{ind_name}-BLAST.report', 'w') as f:
        for peptide, hits in results:
            f.write(f"Peptide: {peptide}\n")
            for hit in hits:
                f.write(f"{hit['hit_id']} - {hit['hit_def']} (E value: {hit['evalue']})\n")
            f.write("=" * 50 + "\n")

def remove_false_AMELY_unique_peptides(ind_name: str, df_fn: str, minimumBLASTLength: int):
    df = pd.read_csv(df_fn)
    
    AMELY_unique_peptides = df[(df['Classification'] == 'AMELY-unique') & (df['Peptide'].str.len() >= minimumBLASTLength)]['Peptide'].tolist()
    if len(AMELY_unique_peptides) > 0:
        print(f"Try to remove AMELY-unique peptides of {ind_name} which have non-AMELY BLAST matching")
        results = autoBLAST(ind_name, AMELY_unique_peptides)
        remove_peptide_list = get_peptides_to_remove(ind_name, results)
        output_BLAST_report(ind_name, results)
        df = df[~df['Peptide'].isin(remove_peptide_list)]
    
    df_out_fn = df_fn.replace('peptides.csv', 'peptides_afterBLAST.csv')
    df.to_csv(df_out_fn, index=False)

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(
        description="Remove AMELY-unique peptides which matching other genes using BLAST",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("ind_name", help="The name of sample")
    parser.add_argument("classified_peptides", help="The path of classified_peptides.csv")
    parser.add_argument("--minimumBLASTLength", "-minL", help="Minimum peptide length for BLAST searches", default=8, type=int)
    args = parser.parse_args()
    
    remove_false_AMELY_unique_peptides(args.ind_name, args.classified_peptides, args.minimumBLASTLength)
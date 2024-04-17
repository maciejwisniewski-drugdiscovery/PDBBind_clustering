import pandas as pd
from pycdhit import cd_hit, read_clstr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os


#  Cluster Output Files
pdbbind_proteins_fasta_output_filepath = '/Users/maciejwisniewski/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta'
cd_hit_output_filepath = '/Users/maciejwisniewski/data/PDBBind_Statistics/Clusters/PDBBind_proteins_cdhit/'

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
dataframe_filepath='/Users/maciejwisniewski/data/PDBBind/LP_PDBBind.csv'
dataframe = pd.read_csv(dataframe_filepath)

# Wczytanie sekwencji białek i SMILES
dataframe['seq'] = dataframe['seq'].str.replace(':','')
smiles = dataframe['smiles'].tolist()
protein_sequences = dataframe['seq'].tolist()
grouped_sequences = {protein_type: dataframe[dataframe['type'] == protein_type]['seq'].tolist()
                     for protein_type in dataframe['type'].unique()}


def merged_fasta(df, protein_merged_seq_output_filepath):

    protein_records = []

    for index, row in df.iterrows():

        pdb_id = row['pdbid']
        description = row['type']

        protein_sequence = row['seq']
        protein_seq_record = SeqRecord(Seq(protein_sequence), id=pdb_id,description=description)
        protein_records.append(protein_seq_record)

    SeqIO.write(protein_records, protein_merged_seq_output_filepath, "fasta")

if not os.path.exists(pdbbind_proteins_fasta_output_filepath):
    merged_fasta(dataframe, pdbbind_proteins_fasta_output_filepath)

res = cd_hit(
    i=pdbbind_proteins_fasta_output_filepath,
    o=cd_hit_output_filepath,
    c=0.7,
    d=0,
    sc=1,
)
a=1
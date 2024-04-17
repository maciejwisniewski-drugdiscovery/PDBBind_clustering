import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import seaborn as sns
from Bio import pairwise2, Align
import numpy as np
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

#  Cluster Output Files
pdbbind_proteins_fasta_output_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta'
cd_hit_output_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_cdhit'
biopython_similarity_matrix_filepath= '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/biopython_pdbbind_sequence_similarity_matrix.npy'

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
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

matrix = Align.substitution_matrices.load("BLOSUM62")

def calculate_similarity(seq1, seq2):
    alignment = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)[0]
    aligned_seq1, aligned_seq2, score, _, _ = alignment
    return score / max(len(seq1), len(seq2))

# Obliczanie macierzy podobieństw
similarity_matrix = np.zeros((len(protein_sequences), len(protein_sequences)))

for i in range(len(protein_sequences)):
    for j in range(i, len(protein_sequences)):
        print(i,':',j)
        similarity_matrix[i][j] = calculate_similarity(protein_sequences[i], protein_sequences[j])
        similarity_matrix[j][i] = similarity_matrix[i][j]


np.save(biopython_similarity_matrix_filepath,similarity_matrix)


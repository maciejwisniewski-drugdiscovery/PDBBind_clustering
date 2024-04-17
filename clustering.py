import os.path

import pandas as pd
from Bio import SeqIO
import pandas as pd
from Bio import pairwise2, Align
import numpy as np

#  Cluster Output Files
pdbbind_proteins_fasta_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta'
clustal_biopython_distance_matrix_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequence_distance_matrix.txt'

cd_hit_output_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_cdhit'

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
dataframe = pd.read_csv(dataframe_filepath)[:3]

# Wczytanie sekwencji białek i SMILES

dataframe['seq'] = dataframe['seq'].str.replace(':','')
smiles = dataframe['smiles'].tolist()
protein_sequences=[]
for index, row in dataframe.iterrows():
    protein_sequences.append(SeqIO.SeqRecord(row['seq'], id=row['pdbid'], description=row['type']))

grouped_sequences = {protein_type: dataframe[dataframe['type'] == protein_type]['seq'].tolist()
                     for protein_type in dataframe['type'].unique()}


# ?
def merged_fasta(df, protein_merged_seq_output_filepath):

    protein_records = []

    for index, row in df.iterrows():

        pdb_id = row['pdbid']
        description = row['type']

        protein_sequence = row['seq']
        protein_seq_record = SeqRecord(Seq(protein_sequence), id=pdb_id,description=description)
        protein_records.append(protein_seq_record)

    SeqIO.write(protein_records, protein_merged_seq_output_filepath, "fasta")


# Zapisujemy sekwencje do pliku w formacie FASTA
if not os.path.exists(pdbbind_proteins_fasta_filepath):
    SeqIO.write(protein_sequences, pdbbind_proteins_fasta_filepath, "fasta")


# Utwórz macierz substytucji BLOSUM62
matrix = Align.substitution_matrices.load("BLOSUM62")

# Funkcja do obliczania podobieństwa między dwiema sekwencjami
def calculate_similarity(seq1, seq2):
    alignment = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5, score_only=True)
    score = alignment
    return score

# Obliczanie macierzy podobieństw
similarity_matrix = np.zeros((len(protein_sequences), len(protein_sequences)))

for i in range(len(protein_sequences)):
    for j in range(i, len(protein_sequences)):
        similarity_matrix[i][j] = calculate_similarity(protein_sequences[i], protein_sequences[j])
        similarity_matrix[j][i] = similarity_matrix[i][j]
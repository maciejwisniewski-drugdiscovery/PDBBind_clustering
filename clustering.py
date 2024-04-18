import os.path

from Bio import SeqIO
import pandas as pd
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq


#  Cluster Output Files
pdbbind_proteins_fasta_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta'
clustalo_distance_matrix_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clustalo_PDBBind_protein_sequences_distance_matrix.txt'
cdhit_clusters_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/cdhit_PDBBind_protein_sequences_clusters.clstr'

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
dataframe = pd.read_csv(dataframe_filepath)
dataframe = dataframe.sort_values(by=['type'])
types = dataframe['type'].drop_duplicates().tolist()
print(types)

# Wczytanie sekwencji białek i SMILES

dataframe['seq'] = dataframe['seq'].str.replace(':','')

smiles = dataframe['smiles'].tolist()
protein_sequences=dataframe['seq'].tolist()

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

    print(protein_records)
    SeqIO.write(protein_records, protein_merged_seq_output_filepath, "fasta")


# Zapisujemy sekwencje do pliku w formacie FASTA
if not os.path.exists(pdbbind_proteins_fasta_filepath):
    merged_fasta(dataframe,pdbbind_proteins_fasta_filepath)



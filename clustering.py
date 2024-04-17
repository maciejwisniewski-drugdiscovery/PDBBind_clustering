import pandas as pd
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO

#  Cluster Output Files
pdbbind_proteins_fasta_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta'
clustal_align_pdbbind_proteins_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_sequences.fasta'

cd_hit_output_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/PDBBind_proteins_cdhit'
biopython_similarity_matrix_filepath= '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/biopython_pdbbind_sequence_similarity_matrix.npy'

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
dataframe = pd.read_csv(dataframe_filepath)

# Wczytanie sekwencji białek i SMILES
dataframe['seq'] = dataframe['seq'].str.replace(':','')
smiles = dataframe['smiles'].tolist()
protein_sequences=[]
for index, row in dataframe.iterrows():
    protein_sequences.append(SeqIO.SeqRecord(row['seq'], id=row['pdbid'], description=row['type']))

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


# Zapisujemy sekwencje do pliku w formacie FASTA
SeqIO.write(protein_sequences, pdbbind_proteins_fasta_filepath, "fasta")

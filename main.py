import pandas as pd
import numpy as np
import os

from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from pycdhit import read_fasta, read_clstr, CDHIT


# Filepath
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/'
cd_hit_directory = '/home2/faculty/mwisniewski/Software/cd-hit-v4.8.1-2019-0228/'


#  Load Raw PDBBind CSV
raw_dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
raw_dataframe = pd.read_csv(raw_dataframe_filepath)




def generate_proteins_fasta(df, proteins_fasta_filepath,pdb_id_column='pdbid',protein_type_column='type',protein_seq_column='seq'):
    '''
    Generate Fasta File for specific proteins from dataframe.
    :param df: Dataframe with PDB ID, Protein_type, Sequence columns
    :param protein_merged_seq_output_filepath:  Output filepath
    :return:
    '''

    protein_records = []

    for index, row in df.iterrows():

        pdb_id = row[pdb_id_column]
        description = row[protein_type_column]
        protein_sequence = row[protein_seq_column]

        protein_seq_record = SeqRecord(Seq(protein_sequence), id=pdb_id,description=description)
        protein_records.append(protein_seq_record)

    SeqIO.write(protein_records, proteins_fasta_filepath, "fasta")
def cluster_proteins_fasta(proteins_fasta_filepath, cd_hit_directory=cd_hit_directory):
    cdhit = CDHIT(prog="cd-hit", path=cd_hit_directory)
    df_in = read_fasta(proteins_fasta_filepath)
    df_out, df_clstr = cdhit.set_options(c=0.95, d=0, n=5).cluster(df_in)
    return df_clstr



#  PIPELINE

# Generate list of protein types
protein_types = raw_dataframe['type'].drop_duplicates().tolist()

# Generate FASTA file for all types
for protein_type in protein_types:
    if not os.path.exists(datadir+'/Clusters/fasta/'+protein_type+'_PDBBind_proteins_sequences.fasta'):
        generate_proteins_fasta(raw_dataframe[raw_dataframe['type'] == protein_type],
                                datadir+'/Clusters/fasta/'+protein_type+'_PDBBind_proteins_sequences.fasta')
        print(protein_type,' - Done')
    else:
        print(protein_type,' - Already exists')

a=1

# CD-HiT clustering
for i, protein_type in enumerate(protein_types):
    protein_type_clusters = cluster_proteins_fasta(datadir+'/Clusters/fasta/'+protein_type+'_PDBBind_proteins_sequences.fasta')

    protein_type_clusters['cluster'] = protein_type_clusters['cluster'].apply(lambda x: protein_type + '_' + str(x))
    protein_type_clusters['identifier'] = protein_type_clusters['identifier'].apply(lambda x: x.split(' ')[0])
    protein_type_clusters.drop(columns=['size', 'identity'],axis=1, inplace=True)

    print(protein_type_clusters)

    if i == 0:
        dataframe = pd.merge(raw_dataframe, protein_type_clusters, left_on='pdbid',right_on='identifier', how='left')
    if i != 0:
        dataframe = pd.merge(dataframe, protein_type_clusters, left_on='pdbid',right_on='identifier', how='left')
        dataframe['cluster_x'] = dataframe['cluster_y'].fillna(dataframe['cluster_x'])
        dataframe['is_representative_x'] = dataframe['is_representative_y'].fillna(dataframe['is_representative_x'])
        dataframe.drop('cluster_y', axis=1, inplace=True)
        dataframe.drop('is_representative_y', axis=1, inplace=True)
        dataframe.drop(columns=['identifier_x'], axis=1, inplace=True)
        dataframe.rename(columns={'cluster_x': 'cluster'}, inplace=True)
        dataframe.rename(columns={'is_representative_x': 'is_representative'}, inplace=True)

    print(dataframe)


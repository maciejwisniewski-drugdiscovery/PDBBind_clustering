import pandas as pd
import numpy as np
#from pycdhit import cd_hit, read_clstr

datadir = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/'

#  Load Raw PDBBind CSV
raw_dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
#raw_dataframe_filepath='/Users/maciejwisniewski/data/PDBBind/LP_PDBBind.csv'

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



#  PIPELINE

# Generate list of protein types
protein_types = raw_dataframe['type'].drop_duplicates().tolist()

# Generate FASTA file for all types
for protein_type in protein_types:
    generate_proteins_fasta(raw_dataframe[raw_dataframe['type'] == protein_type],
                            datadir+'/Clusters/fasta/'+protein_type+'_PDBBind_proteins_sequences.fasta')
a=1

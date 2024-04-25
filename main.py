import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
import torch

from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from pycdhit import read_fasta, CDHIT

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity

# Filepath
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics'
cd_hit_directory = '/home2/faculty/mwisniewski/Software/cd-hit-v4.8.1-2019-0228/'


#  Load Raw PDBBind CSV
raw_dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/ecod_LP_PDBBind.csv'
raw_dataframe = pd.read_csv(raw_dataframe_filepath)
raw_dataframe['seq'] = raw_dataframe['seq'].str.replace(':','')




def generate_proteins_fasta(df, proteins_fasta_filepath,pdb_id_column='pdbid',protein_ecod_column='ECOD_Cluster_4',protein_seq_column='seq'):
    '''
    Generate Fasta File for specific proteins from dataframe.
    :param df: Dataframe with PDB ID, Protein_type, Sequence columns
    :param protein_merged_seq_output_filepath:  Output filepath
    :return:
    '''

    protein_records = []

    for index, row in df.iterrows():

        pdb_id = row[pdb_id_column]
        description = row[protein_ecod_column]
        protein_sequence = row[protein_seq_column]

        protein_seq_record = SeqRecord(Seq(protein_sequence), id=pdb_id,description=description)
        protein_records.append(protein_seq_record)

    SeqIO.write(protein_records, proteins_fasta_filepath, "fasta")
def cluster_proteins_fasta(proteins_fasta_filepath, cd_hit_directory=cd_hit_directory):
    cdhit = CDHIT(prog="cd-hit", path=cd_hit_directory)
    df_in = read_fasta(proteins_fasta_filepath)
    df_out, df_clstr = cdhit.set_options(c=0.95, d=0, n=5).cluster(df_in)
    print(df_clstr)
    return df_clstr
def calculate_SMILES_similarity_matrix(smiles_list):
    # Konwertowanie SMILES na obiekty molekularne
    mols = [Chem.MolFromSmiles(smiles, sanitize=False) for smiles in smiles_list]
    print(mols)

    # Tworzenie fingerprintów Tanimoto
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=7)
    fgrps = [rdkit_gen.GetFingerprint(mol) for mol in mols]

    similarities = np.zeros((len(fgrps), len(fgrps)))
    for i in range(1, len(fgrps)):
        print(i)
        similarity = BulkTanimotoSimilarity(fgrps[i], fgrps[:i])
        similarities[i, :i] = similarity
        similarities[:i, i] = similarity

    return similarities





#  PIPELINE

# Generate list of protein types
protein_types = raw_dataframe['ECOD_Cluster_4'].drop_duplicates().tolist()

# Generate FASTA file for all types
for protein_type in protein_types:
    if not os.path.exists(datadir+'/Clusters/fasta/'+protein_type+'_PDBBind_proteins_sequences.fasta'):
        if len(raw_dataframe[raw_dataframe['ECOD_Cluster_4']==protein_type]) == 1:
            generate_proteins_fasta(raw_dataframe[raw_dataframe['ECOD_Cluster_4'] == protein_type],
                                    datadir+'/Clusters/fasta/'+protein_type.replace('/','').replace('\\','').replace(',','').replace(' ','_').replace('(','').replace(')','')
                                    +'_PDBBind_proteins_sequences.fasta')
            print(protein_type,' - Done')
        else:
            print(protein_type,' - Already exists')


# CD-HiT clustering
if not os.path.exists(datadir+'/Clusters/clusters/cdhit_protein_sequences_clusters.csv'):
    for i, protein_type in enumerate(protein_types):
        print(protein_type)
        if len(raw_dataframe[raw_dataframe['ECOD_Cluster_4'] == protein_type]) == 1:
            protein_type_clusters = cluster_proteins_fasta(datadir+'/Clusters/fasta/'+protein_type.replace('/','').replace('\\','').replace(',','').replace(' ','_').replace('(','').replace(')','')+'_PDBBind_proteins_sequences.fasta')
            protein_type_clusters['cluster'] = protein_type_clusters['cluster'].apply(lambda x: protein_type + '_' + str(x))
            protein_type_clusters['identifier'] = protein_type_clusters['identifier'].apply(lambda x: x.split(' ')[0])
            print(protein_type_clusters)
            protein_type_clusters.drop(columns=['size', 'identity'],axis=1, inplace=True)
        else:
            protein_type_clusters = pd.DataFrame(columns=['identifier', 'cluster','is_representative'])
            temp_row = {'identifier': raw_dataframe[raw_dataframe['ECOD_Cluster_4'] == protein_type]['pdbid'].values[0],
                           'cluster': protein_type,
                           'is_representative': True}
            temp_df = pd.DataFrame(temp_row, index=[0])

            # Dodawanie nowego wiersza do DataFrame
            protein_type_clusters = pd.concat([protein_type_clusters, temp_df], ignore_index=True)
        if i == 0:
            dataframe = pd.merge(raw_dataframe, protein_type_clusters, left_on='pdbid',right_on='identifier', how='left')
        if i != 0:
            dataframe.drop(columns=['identifier'], axis=1, inplace=True)
            dataframe = pd.merge(dataframe, protein_type_clusters, left_on='pdbid',right_on='identifier', how='left')
            dataframe['cluster_x'] = dataframe['cluster_y'].fillna(dataframe['cluster_x'])
            dataframe['is_representative_x'] = dataframe['is_representative_y'].fillna(dataframe['is_representative_x'])
            dataframe.drop('cluster_y', axis=1, inplace=True)
            dataframe.drop('is_representative_y', axis=1, inplace=True)
            dataframe.rename(columns={'cluster_x': 'cluster'}, inplace=True)
            dataframe.rename(columns={'is_representative_x': 'is_representative'}, inplace=True)

    dataframe.drop(columns=['Unnamed: 0.2', 'Unnamed: 0.1', 'Unnamed: 0','identifier'],axis=1,inplace=True)
    print(dataframe.columns)

    # Save Protein Sequence Clusters:
    dataframe = dataframe.sort_values(by=['cluster']).reset_index(drop=True)
    dataframe.to_csv(datadir+'/Clusters/clusters/cdhit_protein_sequences_clusters.csv')
    # Generate Protein Fasta File for ClustalO
    generate_proteins_fasta(dataframe,
                            datadir + '/Clusters/fasta/ClustalO_PDBBind_proteins_sequences.fasta',
                            protein_type_column='cluster')
else:
    dataframe = pd.read_csv(datadir+'/Clusters/clusters/cdhit_protein_sequences_clusters.csv')


# End and run sbatch clustalo_distance_matrix.sh
if not os.path.exists(datadir+'/Clusters/matrices/ClustalO_PDBBind_protein_sequences_distance_matrix.txt'):
    print('Run sbatch clustalo_distance_matrix.sh to continue')
    sys.exit()

# Convert Distance Matrix to Tensor:
if not os.path.exists(datadir+'/Clusters/matrices/ClustalO_PDBBind_protein_sequences_distance_matrix.pt'):
    raw_clustalo_matrix = datadir+'/Clusters/matrices/ClustalO_PDBBind_protein_sequences_distance_matrix.txt'
    numpy_clustalo_matrix = np.genfromtxt(raw_clustalo_matrix, skip_header=1)
    numpy_clustalo_matrix = numpy_clustalo_matrix[:, 1:]
    torch_clustalo_tensor = torch.from_numpy(numpy_clustalo_matrix)
    torch.save(torch_clustalo_tensor, datadir+'/Clusters/matrices/ClustalO_PDBBind_protein_sequences_distance_matrix.pt')
else:
    torch_clustalo_tensor = torch.load(datadir+'/Clusters/matrices/ClustalO_PDBBind_protein_sequences_distance_matrix.pt')

# Generate HeatMap
if not os.path.exists(datadir+'/Clusters/images/ClustalO_PDBBind_protein_sequences_distance_heatmap.png'):

    categories = dataframe['type'].unique()
    divisions = [0] + list(dataframe.groupby('type').size().cumsum())
    print(divisions)

    print('Heat Map Creation')

    plt.figure(figsize=(14, 12))
    plt.imshow(torch_clustalo_tensor, cmap='gist_ncar', interpolation='nearest')

    plt.xticks(np.array(divisions[:-1]) + np.diff(divisions), categories, rotation=45,fontsize=8)
    plt.yticks(np.array(divisions[:-1]) + np.diff(divisions), categories, fontsize=8)

    plt.title('Distance Matrix of Protein Sequences in PDBBind (Clustal Omega)')
    plt.xlabel('Proteins')
    plt.ylabel('Proteins')
    plt.colorbar()

    plt.savefig(datadir+'/Clusters/images/ClustalO_PDBBind_protein_sequences_distance_heatmap.png')

# Smiles Similarity Map
if not os.path.exists(datadir+'/Clusters/matrices/Tanimoto_PDBBind_ligand_SMILES_similarity_matrix.pt'):
    smiles_list = dataframe['smiles'].tolist()
    numpy_tanimoto_matrix = calculate_SMILES_similarity_matrix(smiles_list)
    torch_tanimoto_tensor = torch.from_numpy(numpy_tanimoto_matrix)
    torch.save(torch_tanimoto_tensor, datadir+'/Clusters/matrices/Tanimoto_PDBBind_ligand_SMILES_similarity_matrix.pt')
else:
    torch_tanimoto_tensor = torch.load(datadir+'/Clusters/matrices/Tanimoto_PDBBind_ligand_SMILES_similarity_matrix.pt')

# Generate Tanimoto HeatMap
if not os.path.exists(datadir+'/Clusters/images/Tanimoto_PDBBind_ligand_SMILES_similarity_heatmap.png'):

    categories = dataframe['type'].unique()
    divisions = [0] + list(dataframe.groupby('type').size().cumsum())
    print(divisions)

    print('Heat Map Creation')

    plt.figure(figsize=(14, 12))
    plt.imshow(torch_tanimoto_tensor, cmap='gist_ncar', interpolation='nearest')

    plt.xticks(np.array(divisions[:-1]) + np.diff(divisions), categories, rotation=45,fontsize=8)
    plt.yticks(np.array(divisions[:-1]) + np.diff(divisions), categories, fontsize=8)

    plt.title('Similarity Matrix of Ligand SMILES in PDBBind (Tanimoto FingerPrints)')
    plt.xlabel('Ligands')
    plt.ylabel('Ligands')
    plt.colorbar()

    plt.savefig(datadir+'/Clusters/images/Tanimoto_PDBBind_ligand_SMILES_similarity_heatmap.png')
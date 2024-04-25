import Bio
import numpy as np
import pandas as pd
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path
import os

dataframe_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clusters/cdhit_protein_sequences_clusters.csv'
dataframe = pd.read_csv(dataframe_filepath)
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/lp/protein/pdb'

def calculate_TMScore():
    return None



tm_score_matrix = np.zeros((len(dataframe), len(dataframe)))

for index, rows in dataframe.iterrows():
    protein_1_filepath = os.path.join(datadir,rows['pdbid']+'_protein')
    protein_1 = get_structure(get_pdb_path(protein_1_filepath))

    coords_1, seq_1 = get_residue_data(protein_1)

    for index_2, rows_2 in dataframe[:index+1].iterrows():
        print(index,'--',index_2)

        protein_2_filepath = os.path.join(datadir, rows_2['pdbid'] + '_protein')
        protein_2 = get_structure(get_pdb_path(protein_2_filepath))
        coords_2, seq_2 = get_residue_data(protein_2)

        res = tm_align(coords_1, coords_2, seq_1, seq_2)

        tm_score_matrix[index,index_2] = res.tm_norm_chain1
        tm_score_matrix[index_2,index] = res.tm_norm_chain2

tm_score_tensor = torch.from_numpy(tm_score_matrix)
torch.save(tm_score_tensor, datadir+'/Clusters/matrices/TMScore_PDBBind_protein_structure_score_matrix.pt')




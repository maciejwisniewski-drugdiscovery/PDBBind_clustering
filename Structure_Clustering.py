import Bio
import numpy as np
import pandas as pd
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path
import os

dataframe = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDB_Statistics/Clusters/clusters/cdhit_protein_sequences_clusters.csv'
dataframe = pd.read_csv(dataframe)
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/lp/protein/pdb'

def calculate_TMScore():
    return None



tm_score_matrix = np.zeros((len(dataframe), len(dataframe)))

for index, rows in dataframe.iterrows():
    protein_1_filepath = os.path.join(datadir,rows['pdbid']+'_protein.pdb')
    protein_1 = get_structure(get_pdb_path(protein_1_filepath))
    chain = next(s.get_chains())
    chain = next(s.get_chains())
    for index_2, rows_2 in dataframe[:index+1].iterrows():

        protein_2_filepath = os.path.join(datadir, rows_2['pdbid'] + '_protein.pdb')
        protein_2 = get_structure(get_pdb_path(protein_2_filepath))

        chain = next(s.get_chains())



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import sys

from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

df = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/clusters/cdhit_protein_sequences_clusters.csv'
simi_file='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.csv'
simi_np_file ='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.npy'

df = pd.read_csv(df)

# Similarity Matrix
if not os.path.exists(simi_np_file):
    def smiles_to_fp(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)



    df['Fingerprint'] = df['smiles'].apply(smiles_to_fp)


    similarities = []
    for i in range(len(df)):
        print(i)
        sims=[]
        for j in range(len(df)):
            if df['Fingerprint'][i] is not None and df['Fingerprint'][j] is not None:
                sim = DataStructs.TanimotoSimilarity(df['Fingerprint'][i], df['Fingerprint'][j])
            else:
                sim = None
            sims.append(sim)
        similarities.append(sims)

    print(similarity_df)
    simi_file='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.csv'
    similarity_df.to_csv(simi_file)


    similarity_matrix = similarity_df.to_numpy()
    np.save('/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.npy', similarity_matrix)
else:
    similarity_matrix = np.load(simi_np_file)


# Which out

ligand_similarity_dict = {}
for i, pdb_id in enumerate(df['pdb_id']):
    similar_ligands_indices = np.where(similarity_matrix[i] > 0.95)[0]
    similar_pdb_ids = df_pdb.loc[similar_ligands_indices, 'pdb_id'].tolist()
    ligand_similarity_dict[pdb_id] = similar_pdb_ids
    print(ligand_similarity_dict[pdb_id])



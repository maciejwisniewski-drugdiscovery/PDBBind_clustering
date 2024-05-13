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

df = pd.read_csv(df)

def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

def ligand_clustering():
    return None


df['Fingerprint'] = df['smiles'].apply(smiles_to_fp)
similarities = []
for i in range(len(df)):
    print(i)
    sims = [DataStructs.TanimotoSimilarity(df['Fingerprint'][i], df['Fingerprint'][j]) for j in range(len(df))]
    similarities.append(sims)

print(similarity_df)
simi_file='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.csv'
similarity_df.to_csv(simi_file)


similarity_matrix = similarity_df.to_numpy()
np.save('/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.npy', similarity_matrix)


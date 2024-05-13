import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
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

    similarity_df = pd.DataFrame(similarities, columns=df.index, index=df.index)
    print(similarity_df)
    simi_file='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.csv'
    similarity_df.to_csv(simi_file)


    similarity_matrix = similarity_df.to_numpy()
    np.save('/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/matrices/ligand_similarities.npy', similarity_matrix)
else:
    similarity_df = pd.read_csv(simi_file)
    similarity_matrix = np.load(simi_np_file)

print(similarity_df)
sys.exit()
# Which out

ligand_similarity_dict = {}
ligand_similarity_dict_2 = {}

for i, row in df.iterrows():
    temp_cluster_df = df
    temp_indices = temp_cluster_df.index.tolist()
    temp_matrix = matrix.loc[temp_indices]
    all = temp_matrix[str(i)]
    not_all = all[all>0.95]
    similar_ligands_indices = not_all.index.tolist()
    similar_ligand_ids = df.iloc[similar_ligands_indices]['pdbid'].tolist()
    similar_ligand_ecods = df.iloc[similar_ligands_indices]['ECOD_Cluster_4'].tolist()
    ligand_similarity_dict[row['pdbid']] = {
        similar_ligand_id:similar_ligand_ecod
        for similar_ligand_id,similar_ligand_ecod in zip(similar_ligand_ids,similar_ligand_ecods)
    }
    print(row['pdbid'],":",ligand_similarity_dict[row['pdbid']])


similar_ligands_to_ligand ='/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics/Clusters/groups/ligand_similarities.pkl'
with open(similar_ligands_to_ligand, 'wb') as plik:
    pickle.dump(ligand_similarity_dict, plik)
    # transformer trenowany na decoy'ach i wiedzy z krysztalow ktory jednoczenie przewiduje strukture i binding affinity
    # konsensusowy w tym sensie


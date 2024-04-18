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
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/PDBBind_Statistics'
dataframe = pd.read_csv(datadir+'/Clusters/clusters/cdhit_protein_sequences_clusters.csv')

# Funkcja do obliczania macierzy podobieństwa
def calculate_SMILES_similarity_matrix(smiles_list):
    # Konwertowanie SMILES na obiekty molekularne
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    print(mols)
    # Tworzenie fingerprintów Tanimoto
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
    print(fps)
    # Obliczanie macierzy podobieństwa
    similarity_matrix = GetTanimotoSimMat(fps)
    #similarity_matrix = np.zeros((len(fps), len(fps)))
    #for i in range(len(fps)):
    #    for j in range(len(fps)):
    #        similarity_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    return similarity_matrix


if not os.path.exists(datadir+'/Clusters/matrices/Tanimoto_PDBBind_ligand_SMILES_similarity_matrix.npy'):
    smiles_list = dataframe['smiles'].tolist()
    numpy_tanimoto_matrix = calculate_SMILES_similarity_matrix(smiles_list)
    np.save(datadir+'/Clusters/matrices/Tanimoto_PDBBind_ligand_SMILES_similarity_matrix.npy',numpy_tanimoto_matrix)

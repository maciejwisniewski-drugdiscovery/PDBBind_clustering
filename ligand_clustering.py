import os.path
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np

# Wczytanie DataFrame'u zawierającego SMILES ligandów oraz Sekwencje AA białek kompleksów
#dataframe_filepath='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
dataframe_filepath='/Users/maciejwisniewski/data/PDBBind/LP_PDBBind.csv'
dataframe = pd.read_csv(dataframe_filepath)
dataframe = dataframe.sort_values(by=['type'])
types = dataframe['type'].drop_duplicates().tolist()

# Funkcja do obliczania macierzy podobieństwa
def calculate_similarity_matrix(smiles_list):
    # Konwertowanie SMILES na obiekty molekularne
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    a=1
    # Tworzenie fingerprintów Tanimoto
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]

    # Obliczanie macierzy podobieństwa
    similarity_matrix = np.zeros((len(fps), len(fps)))
    for i in range(len(fps)):
        for j in range(len(fps)):
            similarity_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    return similarity_matrix


similarity_matrix = calculate_similarity_matrix()
print(similarity_matrix)
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser, Structure, NeighborSearch
from collections import Counter
import tempfile
import os
import pandas as pd

dataframe_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
ECOD_dataframe_filepath = ''
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/lp'

dataframe = pd.read_csv(dataframe_filepath)
dataframe['closest_chain_to_ligand'] = ''
dataframe['ECOD'] = ''

def mol2_to_biopython_structure(mol2_file):
    # Wczytanie ligandu z pliku Mol2
    mol = Chem.MolFromMol2File(mol2_file)
    # Tworzenie pliku tymczasowego dla PDB
    temp_pdb_file = tempfile.NamedTemporaryFile(suffix='.pdb',delete=False).name
    # Zapisanie ligandu do pliku PDB
    Chem.MolToPDBFile(mol,temp_pdb_file)
    # Tworzenie parsera PDB
    parser = PDBParser()
    # Parsowanie danych z pliku PDB
    pdb_struct = parser.get_structure('ligand',temp_pdb_file)
    # Usuwanie tymczasowego pliku
    os.unlink(temp_pdb_file)
    return pdb_struct

def find_closest_chain_to_ligand(protein_pdb_file,ligand_mol2_file):
    # Inicjalizacja parsera PDB
    parser = PDBParser()
    # Wczytanie struktur białka i ligandu
    protein_structure = parser.get_structure('protein',protein_pdb_file)
    ligand_structure = mol2_to_biopython_structure(ligand_mol2_file)
    # Pobranie łańcuchów białka
    protein_chains = list(protein_structure.get_chains())
    # Inicjalizacja wyszukiwania sąsiadów
    ns = NeighborSearch(list(protein_structure.get_atoms()))
    ligand_closest_chains = []
    min_distance = float('inf')

    # Iteracja po atomach ligandu
    for ligand_atom in ligand_structure.get_atoms():
        # Znalezienie najbliższego atomu białka dla danego atomu liganda
        closest_atoms = ns.search(ligand_atom.get_coord(),4.8)
        atom_closest_chains = [closest_atom.get_parent() for closest_atom in closest_atoms]
        count_atom_closest_chains = Counter(atom_closest_chains)
        print(count_atom_closest_chains)
        try:
            atom_closest_chain = count_atom_closest_chains.most_common(1)[0][0]
        except:
            atom_closest_chain = None
        ligand_closest_chains.append(atom_closest_chain)

    count_ligand_closest_chains = Counter(ligand_closest_chains)
    ligand_closest_chain = count_ligand_closest_chains.most_common(1)[0][0]
    print(ligand_closest_chain)
    return ligand_closest_chain


dataframe = dataframe[0:1]
for index,row in dataframe.iterrows():

    molecule = row['pdbid']
    print(molecule)
    protein_pdb_file = os.path.join(datadir,'protein','pdb',molecule+'_protein.pdb')
    ligand_mol2_file = os.path.join(datadir,'ligand','mol2',molecule+'_ligand.mol2')
    closest_chain_to_ligand = find_closest_chain_to_ligand(protein_pdb_file,ligand_mol2_file)
    dataframe.at[index,'closest_chain_to_ligand'] = closest_chain_to_ligand

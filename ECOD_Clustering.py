from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser, Structure, NeighborSearch
from collections import Counter
import tempfile
import os
import pandas as pd
pd.set_option('display.max_columns', None)

dataframe_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'
output_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/ecod_LP_PDBBind.csv'

ECOD_dataframe_filepath = '/mnt/evafs/groups/sfglab/mwisniewski/ingenix/data/ECOD/ecod.develop291.domains.txt'
datadir = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/lp'

dataframe = pd.read_csv(dataframe_filepath)
dataframe['ligand_closest_chain'] = ''
dataframe['ligand_closest_residue_id'] = ''
dataframe['ECOD_Cluster_1'] = ''
dataframe['ECOD_Cluster_2'] = ''
dataframe['ECOD_Cluster_3'] = ''
dataframe['ECOD_Cluster_4'] = ''
dataframe['ECOD_Cluster_5'] = ''

ECOD_dataframe = pd.read_csv(ECOD_dataframe_filepath,sep='\t')

def parse_range(s):
    try:
        start, end = map(int, s.split('-'))
        return range(start, end + 1)
    except:
        print(s)

def preprocess_ECOD_df(ECOD_dataframe):
    ECOD_dataframe['pdb_range']
    ECOD_dataframe['chain'] = ECOD_dataframe['chain'].apply(lambda x: str(x).upper())
    ECOD_dataframe = ECOD_dataframe[~ECOD_dataframe['ecod_domain_id'].str.contains('e5j3dA3')]
    ECOD_dataframe['pdb_range'] = ECOD_dataframe['pdb_range'].apply(lambda x: x.split(','))
    ECOD_dataframe = ECOD_dataframe.explode(column=['pdb_range'])
    ECOD_dataframe = ECOD_dataframe[ECOD_dataframe['pdb_range'].str.match(r'^[A-Za-z]+:[0-9]+-[0-9]+$')]
    ECOD_dataframe['pdb_range'] = ECOD_dataframe['pdb_range'].apply(lambda x: x.split(':')[-1])
    ECOD_dataframe['pdb_range'] = ECOD_dataframe['pdb_range'].apply(lambda x: parse_range(x))
    print('\n\n\n')
    print(ECOD_dataframe)
    print('\n\n\n')
    print(len(ECOD_dataframe))
    return ECOD_dataframe


def mol2_to_biopython_structure(mol2_file):
    # Wczytanie ligandu z pliku Mol2
    mol = Chem.MolFromMol2File(mol2_file,sanitize=False)
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
    # Iteracja po atomach ligandu
    for ligand_atom in ligand_structure.get_atoms():
        # Znalezienie najbliższego atomu białka dla danego atomu liganda
        closest_atoms = ns.search(ligand_atom.get_coord(),4.8)
        atom_closest_chains = [(str(closest_atom.get_parent().get_parent())[-2],closest_atom.get_parent().get_id()[1])
                               for closest_atom in closest_atoms
                               if closest_atom.get_parent().get_full_id()[0] == 'protein']
        count_atom_closest_chains = Counter(atom_closest_chains)
        try:
            atom_10_closest_chains = count_atom_closest_chains.most_common(10)
            for atom_closest_chain in atom_10_closest_chains:
                if atom_closest_chain[1] > 4:
                    ligand_closest_chains.append(atom_closest_chain[0])
        except:
            atom_closest_chain = None
    count_ligand_closest_chains = Counter(ligand_closest_chains)
    ligand_closest_chains_and_residues = []
    for ligand_closest_chain_and_residue in count_ligand_closest_chains.most_common(5):
        ligand_closest_chains_and_residues.append(ligand_closest_chain_and_residue[0])


    return ligand_closest_chains_and_residues
def check_range(range_tuple, x):
    try:
        return range_tuple[0] <= x <= range_tuple[-1]
    except:
        False


def find_ECOD(molecule,ligand_closest_chain,ligand_closest_residue_id,ECOD_dataframe):
    option = ECOD_dataframe[ECOD_dataframe['pdb'].str.contains(molecule)]
    if len(option) == 1:
        return (option['arch_name'].values[0],
                option['x_name'].values[0],
                option['h_name'].values[0],
                option['t_name'].values[0],
                option['f_name'].values[0])

    option = option[option['chain'].str.contains(ligand_closest_chain)]
    if len(option) == 1:
        return (option['arch_name'].values[0],
                option['x_name'].values[0],
                option['h_name'].values[0],
                option['t_name'].values[0],
                option['f_name'].values[0])
    option = option[option['pdb_range'].apply(lambda x: len(x) > 0)]
    option = option[option['pdb_range'].apply(lambda r: check_range(r, ligand_closest_residue_id))]
    if len(option) == 1:
        return (option['arch_name'].values[0],
                option['x_name'].values[0],
                option['h_name'].values[0],
                option['t_name'].values[0],
                option['f_name'].values[0])
    else:
        return ('ARCH_UNCLASSIFIED','X_UNCLASSIFIED','H_UNCLASSIFIED','T_UNCLASSIFIED','F_UNCLASSIFIED')

print('ECOD Dataframe Preprocessing')
ECOD_dataframe = preprocess_ECOD_df(ECOD_dataframe)

dataframe = dataframe[dataframe['pdbid']=='1ksn']

for index,row in dataframe.iterrows():
    print(index,'/',len(dataframe))
    molecule = row['pdbid']
    print(molecule)
    protein_pdb_file = os.path.join(datadir,'protein','pdb',molecule+'_protein.pdb')
    ligand_mol2_file = os.path.join(datadir,'ligand','mol2',molecule+'_ligand.mol2')
    ligand_closest_chains_and_residues = find_closest_chain_to_ligand(protein_pdb_file,ligand_mol2_file)
    for ligand_closest_chain_and_residue in ligand_closest_chains_and_residues:
        try:
            ligand_closest_chain,ligand_closest_residue_id = ligand_closest_chain_and_residue
            print('\t',ligand_closest_chain, ligand_closest_residue_id)
            cluster_1, cluster_2, cluster_3, cluster_4, cluster_5 = find_ECOD(molecule, ligand_closest_chain,
                                                                              ligand_closest_residue_id,ECOD_dataframe)
            dataframe.at[index,'ligand_closest_chain'] = ligand_closest_chain
            dataframe.at[index,'ligand_closest_residue_id'] = ligand_closest_residue_id
            print('\t',cluster_1)
            print('\t',cluster_2)
            print('\t',cluster_3)
            print('\t',cluster_4)
            print('\t',cluster_5)
            dataframe.at[index,'ECOD_Cluster_1'] = cluster_1
            dataframe.at[index,'ECOD_Cluster_2'] = cluster_2
            dataframe.at[index,'ECOD_Cluster_3'] = cluster_3
            dataframe.at[index,'ECOD_Cluster_4'] = cluster_4
            dataframe.at[index,'ECOD_Cluster_5'] = cluster_5
            break
        except Exception as e:
            print(e)
#dataframe.to_csv(output_filepath)




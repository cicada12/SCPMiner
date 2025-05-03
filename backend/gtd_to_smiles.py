
from rdkit import Chem
from rdkit.Chem import rdmolops

def graph_transaction_to_smiles(graph_data):
    """
    Converts a single graph transaction into SMILES format.

    Parameters:
    graph_data (str): A string representing the graph transaction for a single molecule.

    Returns:
    str: The SMILES string for the molecule.
    """
    mol = Chem.RWMol()  # RDKit editable molecule object
    atom_indices = {}  # Map of atom indices for the current molecule
    bonds = []  # List to store bond information
    
    for line in graph_data.splitlines():
        line = line.strip()
        
        if line.startswith('v'):
            # Add atom (node)
            _, atom_idx, atom_type = line.split()
            atom_idx = int(atom_idx)
            atom = Chem.Atom(atom_type)
            atom_indices[atom_idx] = mol.AddAtom(atom)  # Add atom to molecule
        
        elif line.startswith('e'):
            # Add bond (edge)
            _, atom1_idx, atom2_idx, bond_type = line.split()
            atom1_idx, atom2_idx = int(atom1_idx), int(atom2_idx)
            bond_type = int(bond_type)  # Convert bond type to integer
            bonds.append((atom1_idx, atom2_idx, bond_type))
    
    # After parsing all atoms and bonds, add bonds to the molecule
    mol = add_bonds_to_molecule(mol, bonds)
    
    # Convert the molecule to SMILES
    smiles = Chem.MolToSmiles(mol)
    return smiles

def add_bonds_to_molecule(mol, bonds):
    """
    Adds bonds to the molecule based on the bond list.
    
    Parameters:
    mol (rdkit.Chem.rdchem.RWMol): The editable RDKit molecule object.
    bonds (list of tuples): A list of tuples where each tuple is (atom1_idx, atom2_idx, bond_type).
    
    Returns:
    mol (rdkit.Chem.rdchem.RWMol): The molecule with added bonds.
    """
    bond_mapping = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE, 4: Chem.BondType.AROMATIC}
    
    for atom1_idx, atom2_idx, bond_type in bonds:
        bond = bond_mapping.get(bond_type, Chem.BondType.SINGLE)  # Default to single bond if unrecognized
        mol.AddBond(atom1_idx, atom2_idx, bond)
    
    return mol

# smiles = graph_transaction_to_smiles(graph_transaction_data)
# print(smiles_to_graph_transaction(smiles_data_act)==smiles_to_graph_transaction(smiles))
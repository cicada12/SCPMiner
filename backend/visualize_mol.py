
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from io import BytesIO

def visualize_molecule(smiles):
    """
    Visualizes the 2D structure of a molecule from a SMILES string and stores the image in a BytesIO buffer.

    Parameters:
    smiles (str): The SMILES string of the molecule.

    Returns:
    img_buffer (BytesIO): A BytesIO buffer containing the image of the molecule.
    """
    # Convert SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles,sanitize=False)
    
    if mol is None:
        print("Invalid SMILES string!")
        return None
    # print()
    # Chem.SanitizeMol(mol)
    # mol = Chem.AddHs(mol)  # Add implicit hydrogens
    AllChem.Compute2DCoords(mol)  # Compute 2D coordinates for visualization
    
    # Create a PIL image object
    img = Draw.MolToImage(mol, size=(720, 720))
    # img.show()
    # Store the image in a BytesIO buffer
    # img_buffer = BytesIO()
    # img.save(img_buffer, format='PNG')  # Save the image to the buffer in PNG format
    # img_buffer.seek(0)  # Reset buffer position to the beginning

    return img

# Example usage
# smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
# img_buffer = visualize_molecule(smiles)

# # If you want to display the image from the buffer, you can use the following code:
# if img_buffer:
#     from PIL import Image
#     img_from_buffer = Image.open(img_buffer)
#     img_from_buffer.show()


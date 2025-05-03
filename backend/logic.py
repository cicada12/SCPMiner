from rdkit import Chem

def run_algorithm(file_path: str, algorithm: str, min_support: float, max_overlap: float, min_coverage: float):
    # Replace this with actual logic
    smiles_to_gspan(file_path)
    print(f"Running algorithm: {algorithm}")
    print(f"File path: {file_path}")
    print(f"Params -> Min Support: {min_support}, Max Overlap: {max_overlap}, Min Coverage: {min_coverage}")
    return {
        "status": "success",
        "details": f"Executed {algorithm} on {file_path}"
    }


def smiles_to_gspan(smiles_list, output_path="dataset.txt"):
    with open(output_path, "w") as f:
        for idx, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue

            f.write(f"t # {idx}\n")
            
            # Atoms
            for atom in mol.GetAtoms():
                f.write(f"v {atom.GetIdx()} {atom.GetSymbol()}\n")
                
            # Bonds
            for bond in mol.GetBonds():
                bgn = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                btype = str(bond.GetBondType())  # e.g. SINGLE, DOUBLE
                f.write(f"e {bgn} {end} {btype}\n")

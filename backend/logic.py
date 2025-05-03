from rdkit import Chem

from gspan import GSpan

def run_algorithm(file_path: str, algorithm: str, min_support: float, max_overlap: float, min_coverage: float):
    with open(file_path, "r") as f:
        smiles_list = [line.strip() for line in f.readlines()]
    smiles_to_gspan(smiles_list)  # creates dataset.txt

    gspan_runner = GSpan("dataset.txt", min_support)
    gspan_runner.mine()

    patterns = gspan_runner.getFrequentSubgraphs()

    return {
        "status": "success",
        "algorithm": "Subgraph Coverage Patterns",
        "parameters": {
            "min_support": min_support,
            "max_overlap": max_overlap,
            "min_coverage": min_coverage
        },
        "executionTime": gspan_runner.getRuntime(),
        "coveredGraphs": gspan_runner.graphCount,
        "discoveredPatterns": patterns
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

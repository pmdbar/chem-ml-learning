import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen

def compute_logp(smiles: str) -> float:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Crippen.MolLogP(mol)

def main():
    df = pd.read_csv("molecules.csv")

    df["logP"] = df["smiles"].apply(compute_logp)

    df.to_csv("molecules_with_logp.csv", index=False)
    print("Wrote molecules_with_logp.csv")

if __name__ == "__main__":
    main()


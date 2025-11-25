import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


def compute_logp(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Crippen.MolLogP(mol)

def compute_mw(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute logP for molecules in a CSV using RDKit."
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Path to input CSV file (must contain a SMILES column).",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Path to output CSV file.",
    )
    parser.add_argument(
        "--smiles-column",
        "-s",
        default="smiles",
        help="Name of the SMILES column in the input CSV (default: smiles).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Reading input from: {args.input}")
    df = pd.read_csv(args.input)

    if args.smiles_column not in df.columns:
        raise ValueError(
            f"SMILES column '{args.smiles_column}' not found in input file. "
            f"Available columns: {list(df.columns)}"
        )
   
    smiles_series = df[args.smiles_column]

    #Compute properties
    df["logP"] = df[args.smiles_column].apply(compute_logp)
    df["mw"] = smiles_series.apply(compute_mw)

    print(f"Writing output to: {args.output}")
    df.to_csv(args.output, index=False)
    print("Done.")


if __name__ == "__main__":
    main()



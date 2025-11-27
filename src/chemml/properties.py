from typing import Optional


from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski


def _mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Parse SMILES into an RDKit Mol or return None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    return mol

def compute_logp(smiles:str) -> Optional[float]:
    """Wildman-Crippen logP, or None if SMILES is invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Crippen.MolLogP(mol)

def compute_mw(smiles: str) -> Optional[float]:
    """Molecular weight (RDKit Descriptors.MolWt), or None if SMILES is invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol)

def compute_canonical_smiles(smiles: str) -> Optional[str]:
    """Canonical SMILES, or None if SMILES is invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)

def compute_tpsa(smiles: str) -> Optional[float]:
    """Topological Polar Surface Area (TPSA), or None if invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Descriptors.TPSA(mol)


#def compute_hbd(smiles: str) -> Optional[int]:
    """Number of H-bond donors, or None if invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Lipinski.NumHDonors(mol)


#def compute_hba(smiles: str) -> Optional[int]:
    """Number of H-bond acceptors, or None if invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Lipinski.NumHAcceptors(mol)


#def compute_rotatable_bonds(smiles: str) -> Optional[int]:
    """Number of rotatable bonds, or None if invalid."""
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Lipinski.NumRotatableBonds(mol)
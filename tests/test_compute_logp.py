import os
import sys
import math

# Add src/ to Python path
current_dir = os.path.dirname(__file__)
project_root = os.path.abspath(os.path.join(current_dir, ".."))
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

from rdkit import Chem
from rdkit.Chem import Crippen
from chemml.compute_logp import compute_logp


def test_compute_logp_matches_rdkit():
    """compute_logp should match RDKit's MolLogP for valid SMILES."""
    smiles = "c1ccccc1"  # benzene
    mol = Chem.MolFromSmiles(smiles)
    rdkit_logp = Crippen.MolLogP(mol)

    my_logp = compute_logp(smiles)

    assert my_logp is not None
    assert math.isclose(my_logp, rdkit_logp, rel_tol=1e-6)


def test_compute_logp_invalid_smiles():
    """Invalid SMILES should return None (no crash)."""
    logp = compute_logp("not_a_smiles")
    assert logp is None


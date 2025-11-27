import math

from rdkit import Chem
from rdkit.Chem import Crippen

from chemml import properties as props  


def test_compute_logp_matches_rdkit():
    smiles = "c1ccccc1"   #benzene
    mol = Chem.MolFromSmiles(smiles)
    rdkit_logp = Crippen.MolLogP(mol)


    my_logp = props.compute_logp(smiles)

    assert my_logp is not None
    assert math.isclose(my_logp, rdkit_logp, rel_tol=1e-6)

def test_compute_logp_invalid_smiles_returns_none():
    logp = props.compute_logp("not_a_smiles")
    assert logp is None


def test_compute_mw_positive_for_valid_smiles():
    mw = props.compute_mw("c1ccccc1")  # benzene
    assert mw is not None
    assert mw > 0


def test_canonical_smiles_roundtrip():
    smiles = "C1=CC=CC=C1"  # non-canonical benzene
    canon = props.compute_canonical_smiles(smiles)
    assert canon == "c1ccccc1"

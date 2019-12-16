"""
Some tests for the compound to reaction search
"""
import pytest
import os
import sys
sys.path.append(os.path.abspath("../../workflow_2"))
import compound_to_reaction as magi
from rdkit import Chem
from rdkit.Chem import AllChem

# Read retro rules database
@pytest.mark.parametrize("diameter, expected_len", [(2,10),(16,5)])
def test_read_retro_rules_db(diameter,expected_len):
    db_path = os.path.join(os.path.dirname(__file__), "./retro_rules_tiny.tsv")
    retro_rules = magi.read_retro_rules(db_path, diameter)
    assert retro_rules.shape == (expected_len,19)

# Mol from SMILES
@pytest.mark.parametrize("smiles", ['Cc1ccccc1', 'Cc1ncccc1'])
def test_good_mol_from_smiles(smiles):
    mol = magi.mol_from_smiles(smiles)
    assert type(mol) == type(Chem.rdchem.Mol())

@pytest.mark.parametrize("smiles", ['not_a_smiles', 23])
def test_bad_mol_from_smiles(smiles):
    try:
        magi.mol_from_smiles(smiles)
        assert False
    except(SystemExit, TypeError):
        assert True

# Test fingerprint similarity
@pytest.mark.parametrize("smile1, smile2, span, expected", 
                        [('Cc1ccccc1', 'Cc1ncccc1', 2, 0.55),
                        ('Cc1ccccc1', 'Cc1ncccc1', 3, 0.5)])
def test_calculate_fingerprint_similarity(smile1, smile2, span, expected):
    mol1 = magi.mol_from_smiles(smile1)
    mol2 = magi.mol_from_smiles(smile2)
    answer = magi.calculate_fingerprint_similarity(mol1, mol2, span)
    assert answer == expected

# Test compound to reaction
def test_compound_to_reaction():
    molecule_smiles = '[CH3][C](=[O])[CH]1[CH2][CH2][CH]2[CH]3[CH2][CH2][CH]4[CH2][C](=[O])[CH2][CH2][C]4([CH3])[CH]3[CH2][CH2][C]12[CH3]'
    db_path = os.path.join(os.path.dirname(__file__), "./retro_rules_c2r_test.tsv")
    rules_to_use = magi.read_retro_rules(db_path, 10)
    fingerprint_radius = 3
    similarity_cutoff = 0.6
    c2r_output_file = "./output.csv"
    expected = "MNXR103630"
    c2r = magi.compound_to_reaction(molecule_smiles, rules_to_use, c2r_output_file,
    fingerprint_radius, similarity_cutoff)
    assert c2r["Reaction_ID"][0] == expected

# Test find matching reactions
@pytest.mark.parametrize("smiles, expected", 
                        [("CCC=O", True),
                        ("CCCO", False)])
def test_mol_matches_reaction(smiles, expected):
    molecule = magi.mol_from_smiles(smiles)
    reaction = AllChem.ReactionFromSmarts("[CH](=[O])-[CH2]-[CH3]>>[CH2](-[OH])-[CH](-[CH3])-[OH]")
    assert magi.mol_matches_reaction(molecule, reaction) == expected

# Test desalting and neutralizing
@pytest.mark.parametrize("smiles, expected", 
[('C(CCN)CC(C(=O)O)N.Cl','NCCCCC(N)C(=O)O'),
('C(CCN)CC(C(=O)O)N','NCCCCC(N)C(=O)O'),
('CC(=O)[O-]','CC(=O)O'),
('CC(=O)O','CC(=O)O'),
('C(CC[NH3+])CC(C(=O)O)N.[Cl-]','NCCCCC(N)C(=O)O'),
('CC(=O)[O-].[Cl-]','CC(=O)O')])
def test_desalt_canonicalize_and_neutralize_smiles(smiles, expected):
    result = magi.desalt_canonicalize_and_neutralize_smiles(smiles)
    assert result == expected

# Test read compounds data 
def test_read_compounds_data():
    df = magi.read_compounds_data(os.path.join(os.path.dirname(__file__), "input_to_prep.csv"))
    assert df.shape == (6,2)

# Test tautomer canonicalization
@pytest.mark.parametrize("smiles, expected", 
[('C/C(O)=N/C','CNC(C)=O'),
('C/C(O)=C/C','CCC(C)=O'),
('N=C1C=CC=CN1','N=c1cccc[nH]1'),
('O/C1=N/CCCCCC1','O=C1CCCCCCN1')])
def test_canonicalize_tautomer(smiles, expected):
    mol = magi.mol_from_smiles(smiles)
    result = magi.canonicalize_tautomer(mol)
    result = Chem.MolToSmiles(result)
    assert result == expected
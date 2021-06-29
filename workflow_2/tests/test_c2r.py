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
    assert len(retro_rules.groups) == (expected_len)

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
                        [('Cc1ccccc1', 'Cc1ncccc1', 2, 0.6944444444444444),
                        ('Cc1ccccc1', 'Cc1ncccc1', 3, 0.5952380952380952)])
def test_calculate_fingerprint_similarity(smile1, smile2, span, expected):
    mol1 = magi.mol_from_smiles(smile1)
    mol2 = magi.mol_from_smiles(smile2)
    answer = magi.calculate_fingerprint_similarity(mol1, mol2, span)
    assert answer == expected

# Test compound to reaction
def test_compound_to_reaction():
    molecule_smiles = '[CH3][C](=[O])[CH]1[CH2][CH2][CH]2[CH]3[CH2][CH2][CH]4[CH2][C](=[O])[CH2][CH2][C]4([CH3])[CH]3[CH2][CH2][C]12[CH3]'
    db_path = os.path.join(os.path.dirname(__file__), "./retro_rules_c2r_test.tsv")
    min_diameter = 10
    rules_to_use = magi.read_retro_rules(db_path, min_diameter)
    fingerprint_radius = 3
    similarity_cutoff = 0.6
    c2r_output_file = "./output.csv"
    expected = "MNXR103630"
    c2r = magi.compound_to_reaction(molecule_smiles = molecule_smiles, 
    rules_to_use = rules_to_use, c2r_output_file = c2r_output_file, min_diameter = min_diameter,
    fingerprint_radius = fingerprint_radius, similarity_cutoff = similarity_cutoff)
    assert c2r["Reaction_ID"][0] == expected

# Test find matching reactions
@pytest.mark.parametrize("smiles, expected", 
                        [("CCC=O", True),
                        ("CCCO", False)])
def test_mol_matches_reaction(smiles, expected):
    molecule = magi.mol_from_smiles(smiles)
    reaction = AllChem.ReactionFromSmarts("[CH](=[O])-[CH2]-[CH3]>>[CH2](-[OH])-[CH](-[CH3])-[OH]")
    assert magi.mol_matches_reaction(molecule, reaction) == expected

@pytest.mark.parametrize("smiles, expected", 
                        [("CCCC(N)CCC", True),
                        ("CCCC(=N)CCC", False)])
def test_mol_matches_reaction2(smiles, expected):
    smiles = magi.prepare_smiles(smiles)
    molecule = magi.mol_from_smiles(smiles)
    reaction = AllChem.ReactionFromSmarts('([#6&v4:1](-[#7&v3:2](-[#1&v1:3])-[#1&v1:4])(-[#6&v4:5])(-[#6&v4:6])-[#1&v1:7])>>([#6&v4:1](-[#6&v4:5])(-[#6&v4:6])=[#8&v2].[#7&v3:2](-[#6&v4](-[#6&v4](-[#6&v4](-[#6&v4](-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1:7])-[#1&v1])(-[#1&v1:3])-[#1&v1:4])')
    assert magi.mol_matches_reaction(molecule, reaction) == expected

# Test desalting and neutralizing
@pytest.mark.parametrize("smiles, expected", 
[('C(CCN)CC(C(=O)O)N.Cl','[H]OC(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]'),
('C(CCN)CC(C(=O)O)N','[H]OC(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]'),
('CC(=O)[O-]','[H]OC(=O)C([H])([H])[H]'),
('CC(=O)O','[H]OC(=O)C([H])([H])[H]'),
('C(CC[NH3+])CC(C(=O)O)N.[Cl-]','[H]OC(=O)C([H])(N([H])[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])N([H])[H]'),
('CC(=O)[O-].[Cl-]','[H]OC(=O)C([H])([H])[H]'),
('[H][O][P](=[O])([O][H])[O][C]1([H])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]1([H])[O][P](=[O])([O][H])[O][H]', 'not_sure')])
def test_prepare_smiles(smiles, expected):
    result = magi.prepare_smiles(smiles)
    assert result == expected

# Test read compounds data 
def test_read_compounds_data():
    df = magi.read_compounds_data(os.path.join(os.path.dirname(__file__), "input_to_prep.csv"))
    assert df.shape == (6,3)

# Test tautomer canonicalization
@pytest.mark.parametrize("smiles, expected", 
[('C/C(O)=N/C','CNC(C)=O'),
('C/C(O)=C/C','CCC(C)=O'),
('N=C1C=CC=CN1','N=c1cccc[nH]1'),
('O/C1=N/CCCCCC1','O=C1CCCCCCN1')])
def test_canonicalize_tautomer(smiles, expected):
    mol1 = Chem.MolFromSmiles(smiles)
    mol2 = Chem.MolFromSmiles(expected)
    result1 = magi.canonicalize_tautomer(mol1)
    result1 = Chem.MolToSmiles(result1)
    result2 = magi.canonicalize_tautomer(mol2)
    result2 = Chem.MolToSmiles(result2)
    assert result1 == result2
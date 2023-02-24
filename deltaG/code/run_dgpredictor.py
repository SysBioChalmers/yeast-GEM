import pandas as pd
import numpy as np
import re
import json
import joblib
import sys

# to your dGPredictor\CC path
sys.path.append(r"D:\model_research\dGPredictor\CC")
sys.path.append('./CC/')
# to your ChemAxon\MarvinSuite\bin\cxcalc.exe path
sys.path.append(r"C:\Program Files\ChemAxon\MarvinSuite\bin\cxcalc.exe")

import chemaxon
from chemaxon import *
from compound import Compound
from compound_cacher import CompoundCacher
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import Draw
from rdkit import Chem

# to your dGPredictor\data\cache_compounds_20160818.csv path
def load_smiles():
    db = pd.read_csv(r'D:\model_research\dGPredictor\data\cache_compounds_20160818.csv',
                     index_col='compound_id')
    db_smiles = db['smiles_pH7'].to_dict()
    return db_smiles

# to your dGPredictor\data\decompose_vector_ac.json path
def load_molsig_rad1():
    molecular_signature_r1 = json.load(open(r'D:\model_research\dGPredictor\data\decompose_vector_ac.json'))
    return molecular_signature_r1

# to your dGPredictor\data\decompose_vector_ac_r2_py3_indent_modified_manual.json path
def load_molsig_rad2():
    molecular_signature_r2 = json.load(
        open(r'D:\model_research\dGPredictor\data\decompose_vector_ac_r2_py3_indent_modified_manual.json'))
    return molecular_signature_r2

# to your dGPredictor\model\M12_model_BR.pkl path
def load_model():
    filename = r'D:\model_research\dGPredictor\model\M12_model_BR.pkl'
    loaded_model = joblib.load(open(filename, 'rb'))
    return loaded_model


def load_compound_cache():
    ccache = CompoundCacher()
    return ccache


def count_substructures(radius, molecule):
    """Helper function for get the information of molecular signature of a
    metabolite. The relaxed signature requires the number of each substructure
    to construct a matrix for each molecule.
    Parameters
    ----------
    radius : int
        the radius is bond-distance that defines how many neighbor atoms should
        be considered in a reaction center.
    molecule : Molecule
        a molecule object create by RDkit (e.g. Chem.MolFromInchi(inchi_code)
        or Chem.MolToSmiles(smiles_code))
    Returns
    -------
    dict
        dictionary of molecular signature for a molecule,
        {smiles: molecular_signature}
    """
    m = molecule
    smi_count = dict()
    atomList = [atom for atom in m.GetAtoms()]

    for i in range(len(atomList)):
        env = Chem.FindAtomEnvironmentOfRadiusN(m, radius, i)
        atoms = set()
        for bidx in env:
            atoms.add(m.GetBondWithIdx(bidx).GetBeginAtomIdx())
            atoms.add(m.GetBondWithIdx(bidx).GetEndAtomIdx())

        # only one atom is in this environment, such as O in H2O
        if len(atoms) == 0:
            atoms = {i}

        smi = Chem.MolFragmentToSmiles(m, atomsToUse=list(atoms),
                                       bondsToUse=env, canonical=True)

        if smi in smi_count:
            smi_count[smi] = smi_count[smi] + 1
        else:
            smi_count[smi] = 1
    return smi_count


def decompse_novel_mets_rad1(novel_smiles, radius=1):
    decompose_vector = dict()

    for cid, smiles_pH7 in novel_smiles.items():
        mol = Chem.MolFromSmiles(smiles_pH7)
        mol = Chem.RemoveHs(mol)
        # Chem.RemoveStereochemistry(mol)
        smi_count = count_substructures(radius, mol)
        decompose_vector[cid] = smi_count
    return decompose_vector


def decompse_novel_mets_rad2(novel_smiles, radius=2):
    decompose_vector = dict()

    for cid, smiles_pH7 in novel_smiles.items():
        mol = Chem.MolFromSmiles(smiles_pH7)
        mol = Chem.RemoveHs(mol)
        # Chem.RemoveStereochemistry(mol)
        smi_count = count_substructures(radius, mol)
        decompose_vector[cid] = smi_count
    return decompose_vector

# def parse_rule(rxn,df_rule):
#     df = df_rule
#     rule_df = df[rxn].to_frame()
#     # new_df = rule_df[(rule_df.T != 0).any()]

#     return rule_df[(rule_df.T != 0).any()]


def parse_reaction_formula_side(s):
    """
        Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        Ignores stoichiometry.

        Returns:
            The set of CIDs.
    """
    if s.strip() == "null":
        return {}

    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if len(tokens) == 0:
            continue
        if len(tokens) == 1:
            amount = 1
            key = member
        else:
            amount = float(tokens[0])
            key = tokens[1]

        compound_bag[key] = compound_bag.get(key, 0) + amount

    return compound_bag


def parse_formula(formula, arrow='<=>', rid=None):
    """
        Parses a two-sided formula such as: 2 C00001 => C00002 + C00003

        Return:
            The set of substrates, products and the direction of the reaction
    """
    tokens = formula.split(arrow)
    if len(tokens) < 2:
        print(('Reaction does not contain the arrow sign (%s): %s'
               % (arrow, formula)))
    if len(tokens) > 2:
        print(('Reaction contains more than one arrow sign (%s): %s'
               % (arrow, formula)))

    left = tokens[0].strip()
    right = tokens[1].strip()

    sparse_reaction = {}
    for cid, count in parse_reaction_formula_side(left).items():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count

    for cid, count in parse_reaction_formula_side(right).items():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count

    return sparse_reaction


def draw_rxn_figure(rxn_dict, db_smiles, novel_smiles):
    # db_smiles = load_smiles()

    left = ''
    right = ''

    for met, stoic in rxn_dict.items():
        if met == "C00080" or met == "C00282":
            continue  # hydogen is not considered
        if stoic > 0:
            if met in db_smiles:
                right = right + db_smiles[met] + '.'
            else:
                right = right + novel_smiles[met] + '.'
        else:
            if met in db_smiles:
                left = left + db_smiles[met] + '.'
            else:
                left = left + novel_smiles[met] + '.'
    smarts = left[:-1] + '>>' + right[:-1]
    # print smarts
    smarts = str(smarts)
    rxn = Reactions.ReactionFromSmarts(smarts, useSmiles=True)
    return Draw.ReactionToImage(rxn)  # , subImgSize=(400, 400))

# def draw_group_changes(rxn,df_rule):
#     df = parse_rule(rxn,df_rule)
#     group_dict = df.to_dict()[rxn]

#     left = ''
#     right = ''

#     for smiles,stoic in group_dict.iteritems():
#         if stoic > 0:
#             right = right + smiles + '.'
#         else:
#             left = left + smiles + '.'
#     smarts = left[:-1] + '>>' + right[:-1]
#     rxn = Reactions.ReactionFromSmarts(smarts, useSmiles=True)
#     return Draw.ReactionToImage(rxn)

# def get_rxn_rule(rid):
#     reaction_dict = json.load(open('../data/optstoic_v3_Sji_dict.json'))
#     molecular_signature = json.load(open('../data/decompose_vector_ac.json'))
#     molsigna_df = pd.DataFrame.from_dict(molecular_signature).fillna(0)
#     all_mets = molsigna_df.columns.tolist()
#     all_mets.append("C00080")
#     all_mets.append("C00282")

#     rule_df = pd.DataFrame(index=molsigna_df.index)

#     info = reaction_dict[rid]

#     # skip the reactions with missing metabolites
#     mets = info.keys()
#     flag = False
#     for met in mets:
#         if met not in all_mets:
#             flag = True
#             break
#     if flag:
#         return None

#     rule_df[rid] = 0
#     for met, stoic in info.items():
#         if met == "C00080" or met == "C00282":
#             continue  # hydogen is zero
#         rule_df[rid] += molsigna_df[met] * stoic
#     return rule_df


def get_rule(rxn_dict, molsig1, molsig2, novel_decomposed1, novel_decomposed2):
    if novel_decomposed1 != None:
        for cid in novel_decomposed1:
            molsig1[cid] = novel_decomposed1[cid]
    if novel_decomposed2 != None:
        for cid in novel_decomposed2:
            molsig2[cid] = novel_decomposed2[cid]

    molsigna_df1 = pd.DataFrame.from_dict(molsig1).fillna(0)
    all_mets1 = molsigna_df1.columns.tolist()
    all_mets1.append("C00080")
    all_mets1.append("C00282")

    molsigna_df2 = pd.DataFrame.from_dict(molsig2).fillna(0)
    all_mets2 = molsigna_df2.columns.tolist()
    all_mets2.append("C00080")
    all_mets2.append("C00282")

# to your path
    moieties_r1 = open(r'D:\model_research\dGPredictor\data\group_names_r1.txt')
    moieties_r2 = open(r'D:\model_research\dGPredictor\data\group_names_r2_py3_modified_manual.txt')
    moie_r1 = moieties_r1.read().splitlines()
    moie_r2 = moieties_r2.read().splitlines()

    molsigna_df1 = molsigna_df1.reindex(moie_r1)
    molsigna_df2 = molsigna_df2.reindex(moie_r2)

    rule_df1 = pd.DataFrame(index=molsigna_df1.index)
    rule_df2 = pd.DataFrame(index=molsigna_df2.index)
    # for rid, value in reaction_dict.items():
    #     # skip the reactions with missing metabolites
    #     mets = value.keys()
    #     flag = False
    #     for met in mets:
    #         if met not in all_mets:
    #             flag = True
    #             break
    #     if flag: continue

    rule_df1['change'] = 0
    for met, stoic in rxn_dict.items():
        if met == "C00080" or met == "C00282":
            continue  # hydogen is zero
        rule_df1['change'] += molsigna_df1[met] * stoic

    rule_df2['change'] = 0
    for met, stoic in rxn_dict.items():
        if met == "C00080" or met == "C00282":
            continue  # hydogen is zero
        rule_df2['change'] += molsigna_df2[met] * stoic

    rule_vec1 = rule_df1.to_numpy().T
    rule_vec2 = rule_df2.to_numpy().T

    m1, n1 = rule_vec1.shape
    m2, n2 = rule_vec2.shape

    zeros1 = np.zeros((m1, 44))
    zeros2 = np.zeros((m2, 44))
    X1 = np.concatenate((rule_vec1, zeros1), 1)
    X2 = np.concatenate((rule_vec2, zeros2), 1)

    rule_comb = np.concatenate((X1, X2), 1)

    # rule_df_final = {}
    # rule_df_final['rad1'] = rule_df1
    # rule_df_final['rad2'] = rule_df2
    return rule_comb, rule_df1, rule_df2


def get_ddG0(rxn_dict, pH, I, novel_mets):
    ccache = CompoundCacher()
    # ddG0 = get_transform_ddG0(rxn_dict, ccache, pH, I, T)
    T = 298.15
    ddG0_forward = 0
    for compound_id, coeff in rxn_dict.items():
        if novel_mets != None and compound_id in novel_mets:
            comp = novel_mets[compound_id]
        else:
            comp = ccache.get_compound(compound_id)
        ddG0_forward += coeff * comp.transform_pH7(pH, I, T)

    return ddG0_forward


def get_dG0(rxn_dict, rid, pH, I, loaded_model, molsig_r1, molsig_r2, novel_decomposed_r1, novel_decomposed_r2, novel_mets):

    # rule_df = get_rxn_rule(rid)
    rule_comb, rule_df1, rule_df2 = get_rule(
        rxn_dict, molsig_r1, molsig_r2, novel_decomposed_r1, novel_decomposed_r2)

    X = rule_comb
    # X = X.reshape(1,-1)
    # pdb.set_trace()
#     print(np.shape(X1))
#     print(np.shape(X2))
#     print(np.shape(X))

    ymean, ystd = loaded_model.predict(X, return_std=True)

    # print(ymean)
    # print(ystd)
    result = {}
    # result['dG0'] = ymean[0] + get_ddG0(rxn_dict, pH, I)
    # result['standard deviation'] = ystd[0]

    # result_df = pd.DataFrame([result])
    # result_df.style.hide_index()
    # return result_df
    return ymean[0] + get_ddG0(rxn_dict, pH, I, novel_mets), ystd[0], rule_df1, rule_df2
    # return ymean[0],ystd[0]


def parse_novel_molecule(add_info):
    result = {}
    for cid, InChI in add_info.items():
        c = Compound.from_inchi('Test', cid, InChI)
        result[cid] = c
    return result


def parse_novel_smiles(result):
    novel_smiles = {}
    for cid, c in result.items():
        smiles = c.smiles_pH7
        novel_smiles[cid] = smiles
    return novel_smiles


if __name__ == '__main__':

    pH = 7.0
    I = 0.1
    db_smiles = load_smiles()
    molsig_r1 = load_molsig_rad1()
    molsig_r2 = load_molsig_rad2()
    loaded_model = load_model()
    # provide the rxn with full kegg annotation
    rxn_str = 'C03044 + C00003 <=> C00810 + C00080 + C00004'
    try:
        rxn_dict = parse_formula(rxn_str)
        mu, std, rule_df1, rule_df2 = get_dG0(
            rxn_dict, 'R00801', pH, I, loaded_model, molsig_r1, molsig_r2, novel_decomposed_r1=None,
            novel_decomposed_r2=None,
            novel_mets=None)
        print('{}: {} kj/mol'.format(rxn_str, mu))
        detaG = mu
    except KeyError as e:
        print('KeyError arised:{} in {}'.format(e, rxn_str))
        detaG = 'KeyError arised:{}'.format(e)
    print(detaG)

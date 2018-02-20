from sympy.utilities.iterables import cartes
import re
from itertools import chain
import numpy as np

def get_elements(formulae):
    elements = []
    for i,f in enumerate(formulae):
        m = re.findall(r'([A-Z][a-z]*)(\d*)', f)
        [elements.append(mm[0]) for mm in m]
    return list(set(elements))

def match_formula_to_elements(r,p):
    elements = get_elements(chain(r,p))
    M = np.zeros((len(elements),sum(1 for _ in chain(p,r))))
    counter = 0
    for f in r:
        m = re.findall(r'([A-Z][a-z]*)(\d*)', f)
        for mm in m:
            try:
                M[elements.index(mm[0]),counter] = float(mm[1]) * -1
            except:
                M[elements.index(mm[0]),counter] = 1.0 * -1
        counter += 1

    for f in p:
        m = re.findall(r'([A-Z][a-z]*)(\d*)', f)
        for mm in m:
            try:
                M[elements.index(mm[0]),counter] = float(mm[1])
            except:
                M[elements.index(mm[0]),counter] = 1.0
        counter += 1        
    return M

def brute_force_balanced_reactions(r,p,max_coefficient=5):
    """
    balance reactions attempting the all vs all combinatorials.
    Setting max_coefficeints of 8 takes ~2 seconds
    Going up to 12 takes about ten seconds

    Args:
        r (iterable): chemical formulae of reactants of a reaction
        p (iterable): chemical formulae of products of a reaction
        max_coefficient (int): maximum allowable coefficient for each compound

    Returns:
        (tuple): coefficients for each compound in order of compounds in r, then p

    Example:
        # for the reaction NH4ClO4 + Al --> Al2O3 + HCl + H2O + N2
        >>> r = ['NH4ClO4', 'Al']
        >>> p = ['Al2O3', 'HCl', 'H2O', 'N2']
        >>> brute_force_balanced_reactions(r, p, max_coefficient=12)
        (6, 10, 5, 6, 9, 3) # (NH4ClO4, Al, Al2O3, HCl, H2O, N2)
    """
    M = match_formula_to_elements(r,p)
#     print(M)
    args = ()
    temp = list(range(1,max_coefficient+1))
    for i in range(M.shape[1]):
        args = args + (temp,)
    counter = 0
    for v in cartes(*args):
        counter += 1
        if np.all(M.dot(v)==0):
            return v
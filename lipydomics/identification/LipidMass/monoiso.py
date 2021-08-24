"""
    LipidMass/monoiso.py
    Dylan H. Ross
    2019/07/22

    description:
        Utilities related to calculating monoisotopic masses
"""


def atomic_mass(atom):
    """
atomic_mass
    description:
        Returns a monoisotopic mass with high precision (6 decimal places) for a specified atom. Reference data taken 
        from https://www.sisweb.com/referenc/source/exactmas.htm
    parameters:
        atom (str) -- atom id 
    returns:
        (float) -- high precision atomic mass 
"""
    # reference dictionary for converting atoms to monoisotopic masses with high precision
    a_to_m = {
        'H': 1.007825,
        'C': 12.,
        'N': 14.003074,
        'O': 15.994915,
        'F': 18.998403,
        'Na': 22.989770,
        'P': 30.973763,
        'S': 31.972072,
        'Cl': 34.968853,
        'K': 38.963708
    }
    if atom not in a_to_m:
        raise ValueError('atomic_mass: atom id {} not available in a_to_m'.format(atom))
    return a_to_m[atom]


def formula_mass(formula):
    """
formula_mass
    description:
        Returns a monoisotopic mass associated with a given chemical formula with high precision (6 decimal places). 
        The formula is provided as a dictionary with atom ids mapped to their counts, i.e. C6H12O6 becomes 
        {'C': 6, 'H': 12, 'O': 6}
    parameters:
        formula (dict(str:int)) -- chemical formula expressed as atom ids mapped to their counts
    returns:
        (float) -- high precision monoisotopic mass
"""
    return sum([atomic_mass(aid) * formula[aid] for aid in formula])


def ms_adduct_mass(neutral_mass, adduct):
    """
Lipid.ms_adduct_mass
    description:
        Returns the monoisotopic m/z with high precision (6 decimal places) of an MS adduct (may be positive or
        negative depending on the net gain in mass from the MS adduct)
    parameters:
        neutral_mass (float) -- neutral mass to combine with adduct mass
        adduct (str) -- MS adduct
    returns:
        (float) -- adduct monoisotopic m/z
"""
    # reference dictionary for converting MS adducts to monoisotopic masses with high precision
    adduct_to_m = {
        '[M]+': {},
        '[M+H]+': {'H': 1},
        '[M+Na]+': {'Na': 1},
        '[M+K]+': {'K': 1},
        '[M+2K]2+': {'K': 2},
        '[M+NH4]+': {'N': 1, 'H': 4},
        '[M+H-H2O]+': {'H': -1, 'O': -1},
        '[M-H]-': {'H': -1},
        '[M+HCOO]-': {'H': 1, 'C': 1, 'O': 2},
        '[M+CH3COO]-': {'H': 3, 'C': 2, 'O': 2},
        '[M-2H]2-': {'H': -2},
        '[M-3H]3-': {'H': -3},
        '[M+Cl]-': {'Cl': 1},
        '[M+2Na-H]+': {'H': -1, 'Na': 2},
        '[M+2H]2+': {'H': 2},
        '[M+3H]3+': {'H': 3}
    }
    if adduct not in adduct_to_m:
        raise ValueError('ms_adduct_mass: MS adduct {} not available in adduct_to_m'.format(adduct))
    # reference dictionary for adduct charges
    adduct_to_z = {
        '[M]+': 1,
        '[M+H]+': 1,
        '[M+Na]+': 1,
        '[M+K]+': 1,
        '[M+2K]2+': 2,
        '[M+NH4]+': 1,
        '[M+H-H2O]+': 1,
        '[M-H]-': -1,
        '[M+HCOO]-': -1,
        '[M+CH3COO]-': -1,
        '[M-2H]2-': -2, 
        '[M-3H]3-': -3,
        '[M+Cl]-': -1,
        '[M+2Na-H]+': 1,
        '[M+2H]2+': 2,
        '[M+3H]3+': 3
    }
    return (formula_mass(adduct_to_m[adduct]) + neutral_mass) / abs(adduct_to_z[adduct])

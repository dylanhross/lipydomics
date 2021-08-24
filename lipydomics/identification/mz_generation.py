"""
    lipydomics/identification/mz_generation.py
    Dylan H. Ross
    2019/10/04

    description:
        utilities to generate predicted m/z data for lipids
"""


from .LipidMass.lipids.glycerolipids import DG, TG
from .LipidMass.lipids.glycolipids import DGDG, GlcADG, MGDG
from .LipidMass.lipids.glycerophospholipids import (
    AcylPG, AcylPE, CL, PA, PC, PE, PG, PI, PIP, PIP2, PIP3, PS, LysylPG, AlanylPG
)
from .LipidMass.lipids.lysoglycerophospholipids import LPA, LPC, LPE, LPG, LPI, LPS, LCL
from .LipidMass.lipids.sphingolipids import Cer, HexCer, GlcCer, SM
from .LipidMass.lipids.sphingolipids.gangliosides import Ganglioside
from .LipidMass.lipids.misc import FA


def enumerate_lipid_class(lipid_class_obj, n_carbon_bounds, n_unsat_bounds, adducts, fa_mod=None, limit_nu=True):
    """
enumerate_lipid_class
    description:
        Generates an exhaustive enumeration of lipid monoisotopic masses corresponding to a single lipid class with
        bounds on the number of carbons and unsaturations as well as a list of MS adducts to cover. The number of
        unsaturations to consider is bounded by carbon chain length:
                  nc < 36 -> max nu = 6
            36 <= nc < 54 -> max nu = 12
            54 <= nc      -> max nu = 18
    parameters:
        lipid_class_obj (LipidMass.Lipid) -- reference to uninitialized lipid class object
        n_carbon_bounds (tuple(int)) -- upper and lower bonds (inclusive) of total carbon number to include
        n_unsat_bounds (tuple(int)) -- upper and lower bonds (inclusive) of total unsaturations to include
        adducts (list(str)) -- all MS adducts to include
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
        [limit_nu (bool)] -- apply the bounds discussed above to the maximum number of unsaturations on the basis of 
                                carbon count [optional, default=True]
    yields:
        (str, str, float) -- name, adduct, monoisotopic mass
"""
    nc_min, nc_max = n_carbon_bounds
    nu_min, nu_max = n_unsat_bounds
    for nc in range(nc_min, nc_max + 1):
        # apply upper limit on unsaturations based on number of carbons
        if limit_nu:
            nu_max_ = min(nu_max, 6) if nc < 36 else (min(nu_max, 12) if nc < 54 else min(nu_max, 18))
        else:
            nu_max_ = nu_max
        for nu in range(nu_min, (nu_max_ + 1)):
            for adduct in adducts:
                lipid = lipid_class_obj(nc, nu, fa_mod=fa_mod) if fa_mod else lipid_class_obj(nc, nu)
                yield lipid.name(), adduct, lipid.ms_adduct_monoiso(adduct)


def enumerate_all_lipids():
    """
enumerate_all_lipids
    description:
        Uses enumerate_lipid_class on all implemented lipid classes to generate all possible lipid masses within a 
        set of constraints on number of carbons and unsaturations
    yields:
        (str, str, float) -- name, adduct, monoisotopic mass
"""
    # general carbon count bounds for diacyl lipids 8,8 to 24,24
    diacyl_nc = (16, 56)
    # general unsaturation bounds for diacyl lipids :0 to :12
    diacyl_nu = (0, 12)

    # diacyl-glyerolipids and glycolipids (DG, MGDG, DGDG, and GlcADG)
    diagls = [DG, MGDG, DGDG, GlcADG]
    diagl_adducts = ['[M+NH4]+', '[M+Na]+', '[M+K]+', '[M-H]-', '[M+HCOO]-', '[M+CH3COO]-', '[M+Cl]-', '[M+H-H2O]+']
    for lc in diagls:
        for l in enumerate_lipid_class(lc, diacyl_nc, diacyl_nu, diagl_adducts):
            yield l

    # MGDG and TG
    for l in enumerate_lipid_class(TG, (24, 72), (0, 24), ['[M+NH4]+', '[M+Na]+', '[M+K]+']):
        yield l

    # diacyl-glycerophospholipids (with plasmalogen and ether derivatives)
    diagpls = [PA, PC, PE, PG, PI, PIP, PIP2, PIP3, PS, LysylPG, AlanylPG]
    diagpl_adducts = ['[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H]-', '[M+HCOO]-', '[M+K]+', '[M+CH3COO]-', '[M+2Na-H]+',
                      '[M+Cl]-']
    for lc in diagpls:
        for fa_mod in [None, 'p', 'o']:
            for l in enumerate_lipid_class(lc, diacyl_nc, diacyl_nu, diagpl_adducts, fa_mod=fa_mod):
                yield l

    # AcylPG, AcylPE, CL, and LCL
    for lc in [AcylPG, AcylPE]:
        for l in enumerate_lipid_class(lc, (24, 64), (0, 18), ['[M-H]-', '[M+Na]+']):
            yield l
    for l in enumerate_lipid_class(CL, (36, 72), (0, 24), ['[M-2H]2-', '[M+2K]2+']):
        yield l
    for l in enumerate_lipid_class(LCL, (24, 64), (0, 18), ['[M-2H]2-', '[M+2K]2+']):
        yield l

    # (monoacyl) lysoglycerophospholipids
    lgpls = [LPA, LPC, LPE, LPG, LPI, LPS]
    lgpl_adducts = ['[M+H]+', '[M+Na]+', '[M+HCOO]-', '[M+CH3COO]-', '[M-H]-']
    for lc in lgpls:
        for fa_mod in [None, 'p', 'o']:
            lgpl_nu = (0, 6)
            for l in enumerate_lipid_class(lc, (12, 24), lgpl_nu, lgpl_adducts, fa_mod=fa_mod):
                yield l

    # sphingolipids (Cer, HexCer, GlcCer, SM)
    sls = [Cer, HexCer, GlcCer, SM]
    sls_adducts = ['[M+H]+', '[M+Na]+', '[M+HCOO]-', '[M+CH3COO]-', '[M-H]-', '[M+K]+', '[M+H-H2O]+']
    for lc in sls:
        # our sphingolipids will all have the 'd' FA mod
        for l in enumerate_lipid_class(lc, (30, 44), (1, 7), sls_adducts, fa_mod='d'):
            yield l

    # gangliosides (only the most common species observed in human brain)
    for gs in ['GA1', 'GA2', 'GA3', 'GM1', 'GM2', 'GM3', 'GD1', 'GD2', 'GD3', 'GT1', 'GT2', 'GT3', 'GQ1']:
        for nc in range(28, 42, 2):
            for nu in range(1, 5):
                for ad in ['[M+H]+', '[M+2H]2+', '[M+3H]3+', '[M-H]-', '[M-2H]2-', '[M-3H]3-']:
                    lipid = Ganglioside(gs, nc, nu, fa_mod='d')
                    yield lipid.name(), ad, lipid.ms_adduct_monoiso(ad)

    # fatty acids
    for l in enumerate_lipid_class(FA, (10, 42), (0, 6), ['[M-H]-']):
        yield l


def get_lipid_mz(lipid_class, lipid_nc, lipid_nu, adduct, fa_mod=None):
    """
get_lipid_mz
    description:
        A helper function for generating a single m/z value from a lipid defined by lipid class, fatty acid sum
        composition, FA modifier, and MS adduct
    parameters:
        lipid_class (str) -- lipid class
        lipid_nc (int) -- sum composition, number of fatty acid carbons
        lipid_nu (int) -- sum composition, number of fatty acid unsaturations
        adduct (str) -- MS adduct
        [fa_mod (None or str)] -- fatty acid modifier (e.g. 'p', 'o') [optional, default=None]
    returns:
        (float) -- predicted m/z
"""
    # map lipid class names to the corresponding generator object
    lipid_classes = {
        # glycerolipids
        'DG': DG, 'TG': TG,
        # glycolipids
        'DGDG': DGDG, 'GlcADG': GlcADG, 'MGDG': MGDG,
        # glycerophospholipids
        'AcylPG': AcylPG, 'AcylPE': AcylPE, 'CL': CL, 'PA': PA, 'PC': PC, 'PE': PE, 'PG': PG, 'PI': PI, 'PIP': PIP,
        'PIP2': PIP2, 'PIP3': PIP3, 'PS': PS, 'LysylPG': LysylPG, 'AlanylPG': AlanylPG,
        # lysoglycerophospholipids
        'LPA': LPA, 'LPC': LPC, 'LPE': LPE, 'LPG': LPG, 'LPI': LPI, 'LPS': LPS, 'LCL': LCL,
        # sphingolipids
        'Cer': Cer, 'HexCer': HexCer, 'GlcCer': GlcCer, 'SM': SM,
        # misc
        'FA': FA
    }

    # make sure the lipid class is defined
    if lipid_class not in lipid_classes:
        m = 'get_lipid_mz: lipid class "{}" not defined'.format(lipid_class)
        raise ValueError(m)

    # predict and return m/z
    lco = lipid_classes[lipid_class]
    lipid = lco(lipid_nc, lipid_nu, fa_mod=fa_mod) if fa_mod else lco(lipid_nc, lipid_nu)
    return lipid.ms_adduct_monoiso(adduct)


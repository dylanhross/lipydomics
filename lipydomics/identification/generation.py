"""
    lipydomics/identification/generation.py
    Dylan H. Ross
    2019/10/04

    description:
        utilities to generate theoretical data for lipids
"""


"""
from LipidMass.lipids.glycerolipids import DG, DGDG, GlcADG, MGDG, TG
from LipidMass.lipids.glycerophospholipids import AcylPG, CL, PA, PC, PE, PG, PI, PS
from LipidMass.lipids.lysoglycerophospholipids import LPA, LPC, LPE, LPG, LPI, LPS
from LipidMass.lipids.sphingolipids import Cer, HexCer, SM
from LipidMass.lipids.misc import FA
"""


def enumerate_lipid_class(lipid_class_obj, n_carbon_bounds, n_unsat_bounds, adducts, fa_mod=None, even_nc_only=False):
    """
enumerate_lipid_class
    description:
        Generates an exhaustive enumeration of lipid monoisotopic masses corresponding to a single lipid class with
        bounds on the number of carbons and unsaturations as well as a list of MS adducts to cover.
    parameters:
        lipid_class_obj (LipidMass.Lipid) -- reference to uninitialized lipid class object
        n_carbon_bounds (tuple(int)) -- upper and lower bonds (inclusive) of total carbon number to include
        n_unsat_bounds (tuple(int)) -- upper and lower bonds (inclusive) of total unsaturations to include
        adducts (list(str)) -- all MS adducts to include
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
    yields:
        (str, str, float) -- name, adduct, monoisotopic mass
"""
    nc_min, nc_max = n_carbon_bounds
    nu_min, nu_max = n_unsat_bounds
    for nc in range(nc_min, nc_max + 1):
        for nu in range(nu_min, nu_max + 1):
            for adduct in adducts:
                lipid = lipid_class_obj(nc, nu, fa_mod=fa_mod) if fa_mod else lipid_class_obj(nc, nu)
                yield lipid.name(), adduct, lipid.ms_adduct_monoiso(adduct)


"""
    LipidMass/base.py
    Dylan H. Ross
    2019/07/22

    description:
        Define base classes
"""


from warnings import warn

from .monoiso import formula_mass, ms_adduct_mass


class Lipid:
    """
Lipid
    description:
        The base class for all lipids.
    attributes:
        lipid_class (str) -- full lipid class (including subclass)
        formula (dict(str:int)) -- chemical formula in the form of a dictionary mapping atom id to atom count
"""

    def __repr__(self):
        """
Lipid.__repr__
    description:
        returns a string representation of this Lipid instance
    returns:
        (str) -- string representation of this Lipid instance
"""
        formula_str = ''
        for a in self.formula:
            formula_str += a
            if self.formula[a] > 1:
                formula_str += str(self.formula[a])
        s = '{}(lipid_class="{}", sum_composition={}, formula="{}")'
        return s.format(self.__class__.__name__,
                        self.lipid_class,
                        self.sum_composition,
                        formula_str)

    def gen_acyl_formula(self, n_acyl, n_carbon, n_unsaturation):
        """
Lipid.gen_acyl_formula
    description:
        Generate a combined chemical formula to account for acyl chains using sum composition
    parameters:
        n_acyl (int) -- number of acyl chains
        n_carbon (int) -- sum FA composition carbons
        n_unsaturation (int) -- sum FA composition unsaturations
    returns:
        (dict(str:int)) -- combined chemical formula for the fatty acyl tails
"""
        c = n_carbon - n_acyl
        h = 2 * c + n_acyl - 2 * n_unsaturation
        return {'C': c, 'H': h}

    def add_to_formula(self, formula):
        """
Lipid.add_to_formula
    description:
        Adds atom counts from a specified chemical formula to this Lipid's chemical formula
    parameters:
        formula (dict(str:int)) -- chemical formula to add
"""
        for a in formula:
            if a not in self.formula:
                self.formula[a] = formula[a]
            else:
                self.formula[a] += formula[a]

    def monoiso(self):
        """
Lipid.monoiso
    description:
        Returns the monoisotopic mass with high precision (6 decimal places)
    returns:
        (float) -- monoisotopic mass
"""
        return round(formula_mass(self.formula), 6)

    def ms_adduct_monoiso(self, adduct):
        """
Lipid.ms_adduct_monoiso
    description:
        Returns the monoisotopic mass with high precision (6 decimal places) of a specified MS adduct
    parameters:
        adduct (str) -- MS adduct
    returns:
        (float) -- monoisotopic mass
"""
        return round(ms_adduct_mass(self.monoiso(), adduct), 6)

    def name(self):
        """
Lipid.name
    description:
        returns a short name for this Lipid instance
    returns:
        (str) -- lipid short name 
"""
        nc, nu = self.sum_composition
        lc = self.lipid_class
        fa_mod = self.fa_mod if hasattr(self, 'fa_mod') else None
        if fa_mod:
            return '{}({}{}:{})'.format(lc, fa_mod, nc, nu)
        else:
            return '{}({}:{})'.format(lc, nc, nu)


class Glycerolipid(Lipid):
    r"""
Glycerolipid
    description:
        Lipid with core structure:

         O __    O--R
        /    \__/
    O==|     /
        \   O
       Ac1   \
              |==O
             /
           Ac2
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
Glycerolipid.__init__
    description:
        Sets the core formula for Glycerolipids and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 5, 'H': 5, 'O': 5}
        # two acyl chains are added to the chemical formula
        acyl_formula = self.gen_acyl_formula(2, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)


class Glycerophospholipid(Lipid):
    r"""
Glycerophospholipid
    description:
        Lipid with core structure:

                   O    R
                   ||  /
         O__    O--P--O
        /   \__/   |
    O==|    /      OH
        \   O
        Ac1  \
              |==O
             /
           Ac2

    * by default the phosphate is protonated, this can be adjusted in head groups with permanent charges by
      subtracting a single H in order to make the lipid neutral *

    * ether or plasmalogen derivatives may be specified with the 'o' or 'p' fa_mod, respectively *
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
Glycerophospholipid.__init__
    description:
        Sets the core formula for Glycerophospholipids and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 5, 'H': 6, 'O': 8, 'P': 1}
        # two acyl chains are added to the chemical formula
        acyl_formula = self.gen_acyl_formula(2, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)
        self.fa_mod = fa_mod
        if fa_mod:
            if fa_mod == 'o':
                # remove 1 O and add 2 H, otherwise behavior is the same
                self.add_to_formula({'O': -1, 'H': 2})
            elif fa_mod == 'p':
                # remove 1 oxygen, make sure there is at least 1 unsaturation
                self.add_to_formula({'O': -1})
                # no more warning for 0 unsaturations in plasmalogen, LIPID MAPS has them
                #if sum_unsaturation < 1:
                #    msg = 'Glycerophospholipid: __init__: sum composition ({}:{}) is unusual for plasmalogen'
                #    warn(msg.format(sum_carbon, sum_unsaturation))


class Lysoglycerophospholipid(Lipid):
    r"""
Lysoglycerophospholipid
    description:
        Lipid with core structure:

                   O    R
                   ||  /
         O__    O--P--O
        /   \__/   |
    O==|    /      OH
        \   OH
        Ac1

    * by default the phosphate is protonated, this can be adjusted in head groups with permanent charges by
      subtracting a single H in order to make the lipid neutral *
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
Lysoglycerophospholipid.__init__
    description:
        Sets the core formula for Lysoglycerophospholipids and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 4, 'H': 7, 'O': 7, 'P': 1}
        # two acyl chains are added to the chemical formula
        acyl_formula = self.gen_acyl_formula(1, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)
        self.fa_mod = fa_mod
        if fa_mod:
            if fa_mod == 'o':
                # remove 1 O and add 2 H, otherwise behavior is the same
                self.add_to_formula({'O': -1, 'H': 2})
            elif fa_mod == 'p':
                # remove 1 oxygen, make sure there is at least 1 unsaturation
                self.add_to_formula({'O': -1})
                # no more warning for 0 unsaturations in plasmalogen, LIPID MAPS has them
                #if sum_unsaturation < 1:
                #    msg = 'Lysoglycerophospholipid: __init__: sum composition ({}:{}) is unusual for plasmalogen'
                #    warn(msg.format(sum_carbon, sum_unsaturation))


class Sphingolipid(Lipid):
    r"""
Sphingolipid
    description:
        lipids with core structure:

                    Ac1
                     \
        C12__         |==O
             \       /
              ==    NH
                \__/
                /  \__O--R
              HO
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
Sphingolipid.__init__
    description:
        Sets the core formula for Sphingolipids and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
        [fa_mod (None or str)] -- fatty acid modifier or None [optional, default=None]
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 19, 'H': 35, 'O': 3, 'N': 1}
        # only one acyl chain is added to the chemical formula, AND 18 carbons and 1 unsaturation are already included
        # in the sphingosine backbone, so these values are subtracted from the sum composition
        sum_carbon_new = sum_carbon - 18
        sum_unsaturation_new = sum_unsaturation - 1
        if sum_carbon_new < 1 or sum_unsaturation_new < 0:
            msg = 'Sphingolipid: __init__: sum composition ({}:{}) is unusual for sphingolipids with sphingosine ' + \
                  'backbone (18:1) and may not be valid'
            warn(msg.format(sum_carbon, sum_unsaturation))

        acyl_formula = self.gen_acyl_formula(1, sum_carbon_new, sum_unsaturation_new)
        self.add_to_formula(acyl_formula)
        self.fa_mod = fa_mod
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)

"""
    LipidMass/lipids/misc.py
    Dylan H. Ross
    2019/07/23

    description:
        Define miscellaneous lipids that do not fit into the other superclasses
"""


from ..base import Lipid


class FA(Lipid):
    r"""
FA
    description:
        free fatty acid:

              O
        Ac1__//
             \
              OH
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
FA.__init__
    description:
        Sets the core formula for FAs and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
"""
        # formula is initially set to the core formula (only carboxylic acid)
        self.formula = {'C': 1, 'H': 1, 'O': 2}
        # one acyl chain is added to the chemical formula
        acyl_formula = self.gen_acyl_formula(1, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)
        self.lipid_class = 'FA'

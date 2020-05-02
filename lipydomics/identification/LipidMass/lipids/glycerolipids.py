"""
    LipidMass/lipids/glycerolipids.py
    Dylan H. Ross
    2019/07/23

    description:
        Define individual glycerolipid classes
"""

from ..base import Glycerolipid


class DG(Glycerolipid):
    """
DG
    description:
        Diacylglycerol lipid class
        Glycerolipid with R-group =

          R  =  --H
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
DG.__init__
    description:
        Initializes an instance of a DG lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation)
        # the R-group for a DG is an --H
        self.add_to_formula({'H': 1})
        # the lipid class is DG
        self.lipid_class = 'DG'


class TG(Glycerolipid):
    r"""
TG
    description:
        Triacylglycerol
        Glycerolipid with core structure:

                 Ac3
                  \
                   |==O
                  /
         O __    O
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
TG.__init__
    description:
        Sets the core formula for TGs and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 6, 'H': 5, 'O': 6}
        # two acyl chains are added to the chemical formula
        acyl_formula = self.gen_acyl_formula(3, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)
        self.lipid_class = 'TG'


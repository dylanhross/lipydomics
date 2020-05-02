"""
    LipidMass/lipids/glycolipids.py
    Dylan H. Ross
    2020/05/02

    description:
        Define individual glycolipid classes
"""

from ..base import Glycerolipid


class MGDG(Glycerolipid):
    r"""
MGDG
    description:
        Monogalactosyl-diacylglycerol lipid class
        Glycolipid with R-group =

               HO      OH
                 \____/
                 /    \
        R  =  __/      \__OH
                \      /
                 O____/
                      \__OH
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
MGDG.__init__
    description:
        Initializes an instance of an MGDG lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation)
        # the R-group for a MGDG is a galactose (C6H11O5)
        self.add_to_formula({'C': 6, 'H': 11, 'O': 5})
        # the lipid class is MGDG
        self.lipid_class = 'MGDG'


class GlcADG(Glycerolipid):
    r"""
GlcADG
    description:
        Glucuronosyl-diacylglycerol lipid class
        Glycolipid with R-group =

               HO      OH
                 \____/
                 /    \
        R  =  __/      \__OH
                \      /
                 O____/
                      \
                       |==O
                      /
                    HO
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
GlcADG.__init__
    description:
        Initializes an instance of an GlcADG lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation)
        # the R-group for a GlcADG is a glucuronic acid (C6H9O6)
        self.add_to_formula({'C': 6, 'H': 9, 'O': 6})
        # the lipid class is GlcADG
        self.lipid_class = 'GlcADG'


class DGDG(Glycerolipid):
    r"""
DGDG
    description:
        Digalactosyl-diacylglycerol lipid class
        Glycolipid with R-group =

               HO      OH
                 \____/
                 /    \
        R  =  __/      \__OH  OH    OH
                \      /      \____/
                 O____/       /    \
                      \__O___/      \__OH
                             \      /
                              O____/
                                   \__OH
 """

    def __init__(self, sum_carbon, sum_unsaturation):
        """
DGDG.__init__
    description:
        Initializes an instance of a DGDG lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation)
        # the R-group for a DGDG is a digalactose (C12H21O10)
        self.add_to_formula({'C': 12, 'H': 21, 'O': 10})
        # the lipid class is DGDG
        self.lipid_class = 'DGDG'
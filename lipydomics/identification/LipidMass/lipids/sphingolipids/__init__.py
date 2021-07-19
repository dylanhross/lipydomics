"""
    LipidMass/lipids/sphingolipids/__init__.py
    Dylan H. Ross
    2019/07/23

    description:
        Define individual sphingolipid classes
"""


from ...base import Sphingolipid


class Cer(Sphingolipid):
    """
Cer
    description:
        Ceramide lipid class
        Sphingolipid with R-group =

          R  =  --H
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
Cer.__init__
    description:
        Initializes an instance of a Cer lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a Cer is an --H
        self.add_to_formula({'H': 1})
        # the lipid class is Cer
        self.lipid_class = 'Cer'


class SM(Sphingolipid):
    r"""
SM
    description:
        Sphingomyelin lipid class
        Sphingolipid with R-group =

                 O(-)
                 |
         R  =  --P==O
                 |
                 O__    |
                    \__N(+)__
                        |
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
SM.__init__
    description:
        Initializes an instance of a SM lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a SM is a phosphocholine group
        self.add_to_formula({'C': 5, 'H': 13, 'O': 3, 'N': 1, 'P': 1})
        # the lipid class is SM
        self.lipid_class = 'SM'


class HexCer(Sphingolipid):
    r"""
HexCer
    description:
        Hexosyl-Ceramide lipid class
        Sphingolipid with R-group =

               HO      OH
                 \____/
                 /    \
        R  =  __/      \__OH
                \      /
                 O____/
                      \__OH
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
HexCer.__init__
    description:
        Initializes an instance of a HexCer lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a HexCer is a hexose (C6H11O5)
        self.add_to_formula({'C': 6, 'H': 11, 'O': 5})
        # the lipid class is HexCer
        self.lipid_class = 'HexCer'


class GlcCer(Sphingolipid):
    r"""
GlcCer
    description:
        Glucosyl-Ceramide lipid class
        Sphingolipid with R-group =

               HO      OH
                 \____/
                 /    \
        R  =  __/      \__OH
                \      /
                 O____/
                      \__OH
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
GlcCer.__init__
    description:
        Initializes an instance of a GlcCer lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a HexCer is a hexose (C6H11O5)
        self.add_to_formula({'C': 6, 'H': 11, 'O': 5})
        # the lipid class is HexCer
        self.lipid_class = 'GlcCer'

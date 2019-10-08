"""
    LipidMass/lipids/lysoglycerophospholipids.py
    Dylan H. Ross
    2019/07/23

    description:
        Define individual lysoglycerophospholipid classes
"""


from LipidMass.base import Lysoglycerophospholipid


class LPA(Lysoglycerophospholipid):
    """
LPA
    description:
        Lysophosphatidic acid lipid class
        Lysoglycerophospholipid with R-group =

         R  =  --H
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
LPA.__init__
    description:
        Initializes an instance of a LPA lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a LPA is an --H
        self.add_to_formula({'H': 1})
        # the lipid class is LPA
        self.lipid_class = 'LPA'


class LPC(Lysoglycerophospholipid):
    r"""
LPC
    description:
        Lysophosphatidylcholine lipid class
        Lysoglycerophospholipid with R-group =

        R  =  __    |
                \__N(+)__
                    |
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
LPC.__init__
    description:
        Initializes an instance of a LPC lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a LPC is a choline, because of the permanent charge on the choline group, 1 H must be
        # subtracted from the formula to account for a negative charge on the phosphate group to make the lipid
        # have an overall neutral charge
        self.add_to_formula({'C': 5, 'H': 12, 'N': 1})
        # the lipid class is LPC
        self.lipid_class = 'LPC'


class LPE(Lysoglycerophospholipid):
    r"""
LPE
    description:
        Lysophosphatidylethanolamine lipid class
        Lysoglycerophospholipid with R-group =

        R  =  __
                \__NH2

"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
LPE.__init__
    description:
        Initializes an instance of a LPE lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a LPE is an ethanolamine (C2H6N)
        self.add_to_formula({'C': 2, 'H': 6, 'N': 1})
        # the lipid class is LPE
        self.lipid_class = 'LPE'


class LPG(Lysoglycerophospholipid):
    r"""
LPG
    description:
        Lysophosphatidylglycerol lipid class
        Lysoglycerophospholipid with R-group =

        R  =  __    OH
                \__/
                /
              HO
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
LPG.__init__
    description:
        Initializes an instance of a LPG lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a LPG is an glycerol (C3H7O2)
        self.add_to_formula({'C': 3, 'H': 7, 'O': 2})
        # the lipid class is LPG
        self.lipid_class = 'LPG'


class LPI(Lysoglycerophospholipid):
    r"""
LPI
    description:
        Lysophosphatidylinositol lipid class
        Lysoglycerophospholipid with R-group =

               HO     OH
                 \___/
                 /   \
        R  =  __/     \__OH
                \     /
                 \___/
                 /   \
               HO     OH
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
LPI.__init__
    description:
        Initializes an instance of a LPI lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a LPI is an inositol (C6H11O5)
        self.add_to_formula({'C': 6, 'H': 11, 'O': 5})
        # the lipid class is LPI
        self.lipid_class = 'LPI'


class LPS(Lysoglycerophospholipid):
    r"""
LPS
    description:
        Lysophosphatidylserine lipid class
        Lysoglycerophospholipid with R-group =

        R  =  __
                \__NH2
                /
            O==|
                \
                 OH
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
LPS.__init__
    description:
        Initializes an instance of a LPS lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a LPS is a serine (C3H6O2N)
        self.add_to_formula({'C': 3, 'H': 6, 'O': 2, 'N': 1})
        # the lipid class is LPS
        self.lipid_class = 'LPS'


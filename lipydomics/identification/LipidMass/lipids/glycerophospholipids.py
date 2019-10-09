"""
    LipidMass/lipids/glycerophospholipids.py
    Dylan H. Ross
    2019/07/23

    description:
        Define individual glycerophospholipid classes
"""


from LipidMass.base import Glycerophospholipid


class PA(Glycerophospholipid):
    """
PA
    description:
        phosphatidic acid lipid class
        Glycerophospholipid with R-group =

         R  =  --H
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
PA.__init__
    description:
        Initializes an instance of a PA lipid
    parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a PA is an --H
        self.add_to_formula({'H': 1})
        # the lipid class is PA
        self.lipid_class = 'PA'


class PC(Glycerophospholipid):
    r"""
PC
    description:
        phosphatidylcholine lipid class
        Glycerophospholipid with R-group =

        R  =  __    |
                \__N(+)__
                    |
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
PC.__init__
    description:
        Initializes an instance of a PC lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a PC is a choline, because of the permanent charge on the choline group, 1 H must be
        # subtracted from the formula to account for a negative charge on the phosphate group to make the lipid
        # have an overall neutral charge
        self.add_to_formula({'C': 5, 'H': 12, 'N': 1})
        # the lipid class is PC
        self.lipid_class = 'PC'


class PE(Glycerophospholipid):
    r"""
PE
    description:
        phosphatidylethanolamine lipid class
        Glycerophospholipid with R-group =

        R  =  __
                \__NH2

"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
PE.__init__
    description:
        Initializes an instance of a PE lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a PE is an ethanolamine (C2H6N)
        self.add_to_formula({'C': 2, 'H': 6, 'N': 1})
        # the lipid class is PE
        self.lipid_class = 'PE'


class PG(Glycerophospholipid):
    r"""
PG
    description:
        phosphatidylglycerol lipid class
        Glycerophospholipid with R-group =

        R  =  __    OH
                \__/
                /
              HO
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
PG.__init__
    description:
        Initializes an instance of a PG lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a PG is an glycerol (C3H7O2)
        self.add_to_formula({'C': 3, 'H': 7, 'O': 2})
        # the lipid class is PG
        self.lipid_class = 'PG'


class PI(Glycerophospholipid):
    r"""
PI
    description:
        phosphatidylinositol lipid class
        Glycerophospholipid with R-group =

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
PI.__init__
    description:
        Initializes an instance of a PI lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a PI is an inositol (C6H11O5)
        self.add_to_formula({'C': 6, 'H': 11, 'O': 5})
        # the lipid class is PI
        self.lipid_class = 'PI'


class PS(Glycerophospholipid):
    r"""
PS
    description:
        Phosphatidylserine lipid class
        Glycerophospholipid with R-group =

        R  =  __
                \__NH2
                /
            O==|
                \
                 OH
"""

    def __init__(self, sum_carbon, sum_unsaturation, fa_mod=None):
        """
PS.__init__
    description:
        Initializes an instance of a PS lipid
    Parameters:
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier to indicate plasmalogen or ether lipids ('p' and 'o', 
                                    respectively) or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group for a PS is a serine (C3H6O2N)
        self.add_to_formula({'C': 3, 'H': 6, 'O': 2, 'N': 1})
        # the lipid class is PS
        self.lipid_class = 'PS'


class CL(Glycerophospholipid):
    r"""
CL
    description:
        Cardiolipin
        Glycerophospholipid with core structure:

                                OH
                                |
                  O   O__    O--P==O
                  \\ /   \__/   |
         O__    O--P     /      O    __O
        /   \__/   |   HO        \__/   \
    O==|    /      OH               \    |==O
        \   O                        O  /
        Ac1  \                      /  Ac4
              |==O              O==|
             /                      \
           Ac2                      Ac3
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
CL.__init__
    description:
        Sets the core formula for Cardiolipins and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 13, 'H': 18, 'O': 17, 'P': 2}
        # two acyl chains are added to the chemical formula
        acyl_formula = self.gen_acyl_formula(4, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)
        self.lipid_class = 'CL'


class AcylPG(Glycerophospholipid):
    r"""
acylPG
    description:
        acylphosphatidylglycerol
        Glycerophospholipid with core structure:


                                   
                   O    
                   ||            O
         O__    O--P--O__    O__//    
        /   \__/   |     \__/   \
    O==|    /      OH    /      Ac3
        \   O          HO
        Ac1  \
              |==O
             /
           Ac2
"""

    def __init__(self, sum_carbon, sum_unsaturation):
        """
acylPG.__init__
    description:
        Sets the core formula for Cardiolipins and determines the formula contribution from the acyl chains
    parameters:
        sum_carbon (int) -- sum acyl carbons
        sum_unsaturation (int) -- sum acyl unsaturations
"""
        # formula is initially set to the core formula (only partial head group, no fatty acid tails)
        self.formula = {'C': 9, 'H': 12, 'O': 11, 'P': 1}
        # two acyl chains are added to the chemical formula
        acyl_formula = self.gen_acyl_formula(3, sum_carbon, sum_unsaturation)
        self.add_to_formula(acyl_formula)
        # store the sum composition
        self.sum_composition = (sum_carbon, sum_unsaturation)
        self.lipid_class = 'AcylPG'

"""
    LipidMass/lipids/sphingolipids/gangliosides.py
    Dylan H. Ross
    2021/07/16

    description:
        Define ganglioside lipid classes
"""


from ...base import Sphingolipid


class Ganglioside(Sphingolipid):
    r"""
Ganglioside
    description:
        lipid class covering gangliosides. Uses standard ganglioside nomenclature (Gxy[z]) where x denotes how many
        sialic acids are attached (A = 0, M = 1, D = 2, T = 3, Q = 4, P = 5), y denotes the nuclear core structure 
        (1 = Galβ1-3GalNAcβ1-4Galβ1-4Glcβ1-1’Cer, 2 = lacking terminal Gal, 3 = lacking terminal Galβ1-3GalNAc, 
        4 = Galβ1-3GalNAcβ1-4Gal) and z denotes how many sialic acids are attached at the interior Gal which is only 
        applicable for nuclear core structure 1 with terminal and internal Gal and is ignored in this case because it
        will not change the results.

    ref:
        Sipione S, Monyror J, Galleguillos D, Steinberg N and Kadam V (2020) Gangliosides in the Brain: Physiology,
        Pathophysiology and Therapeutic Applications. Front. Neurosci. 14:572965. doi: 10.3389/fnins.2020.572965
"""

    @staticmethod
    def ganglioside_to_r_group(ganglioside):
        """
Ganglioside.ganglioside_to_r_group
    description:
        returns the R-group formula corresponding to a ganglioside name
    parameters:
        ganglioside (str) -- specify the ganglioside with standard naming convention (Gxy[z])
    returns:
        (dict(str:int)) -- chemical formula of R-group
"""
        # parse x and y from the ganglioside name, make sure they are valid
        x, y = ganglioside[1], ganglioside[2]
        if y not in ['1', '2', '3', '4']:
            msg = 'Ganglioside: ganglioside_to_r_group: ganglioside ' \
                  '"{}" not recognized, y should be 1-4'.format(ganglioside)
            raise ValueError(msg)
        if x not in ['A', 'M', 'D', 'T', 'Q', 'P', 'H', 'S', 'O']:
            msg = 'Ganglioside: ganglioside_to_r_group: ganglioside ' \
                  '"{}" not recognized, x should be A, M, D, T, Q, P, H, S or O'.format(ganglioside)
            raise ValueError(msg)
        # use y first to determine the nuclear core structure
        r = {'C': 6, 'H': 11, 'O': 5, 'N': 0}
        i = 4 - int(y)
        while i > 0:
            if i == 2:
                r['C'] += 8
                r['H'] += 13
                r['O'] += 5
                r['N'] += 1
            else:
                r['C'] += 6
                r['H'] += 10
                r['O'] += 5
            i -= 1
        # add sialic acids
        for _ in range({'A': 0, 'M': 1, 'D': 2, 'T': 3, 'Q': 4, 'P': 5, 'H': 6, 'S': 7, 'O': 8}[x]):
            r['C'] += 11
            r['H'] += 17
            r['O'] += 8
            r['N'] += 1
        # get rid of N in the formula if none were added
        if r['N'] == 0:
            _ = r.pop('N', None)
        return r

    def __init__(self, ganglioside, sum_carbon, sum_unsaturation, fa_mod=None):
        """
Ganglioside.__init__
    description:
        Initializes an instance of a GlcCer lipid
    parameters:
        ganglioside (str) -- specify the ganglioside with standard naming convention (Gxy[z])
        sum_carbon (int) -- sum FA composition carbons
        sum_unsaturation (int) -- sum FA composition unsaturations
        [fa_mod (None or str)] -- fatty acid modifier or None [optional, default=None]
"""
        super().__init__(sum_carbon, sum_unsaturation, fa_mod=fa_mod)
        # the R-group is determined by the ganglioside name
        self.add_to_formula(self.ganglioside_to_r_group(ganglioside))
        # lipid class is the ganglioside type
        self.lipid_class = ganglioside

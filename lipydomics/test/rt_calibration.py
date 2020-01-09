"""
    lipydomics/test/rt_calibration.py
    Dylan H. Ross
    2020/01/09

    description:
        Define tests associated with the lipydomics/identification/rt_calibration module
"""


import os

from lipydomics.identification.rt_calibration import get_ref_rt
from lipydomics.identification.lipid_parser import parse_lipid


def get_ref_rt_lipids1():
    """
get_ref_rt_real1
    description:
        Tests the get_ref_rt helper function using the lipids from lipids_1.txt using both strict and non-strict modes

        Test fails if there are any errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    with open(os.path.join(os.path.dirname(__file__), 'lipids_1.txt'), 'r') as f:
        for line in f.readlines():
            parsed = parse_lipid(line.strip())
            fam = parsed['fa_mod'] if 'fa_mod' in parsed else None
            rt_strict = get_ref_rt(parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat'], fa_mod=fam)
            rt = get_ref_rt(parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat'], fa_mod=fam, strict=False)

    return True


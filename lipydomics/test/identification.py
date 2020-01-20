"""
    lipydomics/test/identification.py
    Dylan H. Ross
    2019/11/26

    description:
        Define tests associated with the lipydomics/identification module
"""


import os

from lipydomics.data import Dataset
from lipydomics.identification import add_feature_ids


def add_feature_ids_any_real1():
    """
add_feature_ids_any_real1
    description:
        Uses the raw data from real_data_1.csv to make compound identifications at any level of confidence. This should
        be a good all-around test for the various identification functions since there are several features in this
        dataset that should not be able to be identified at any level (and thus will run through all of the functions).

        Test fails if there are any errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')

    add_feature_ids(dset, [0.05, 0.5, 5.])

    return True

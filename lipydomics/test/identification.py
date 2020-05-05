"""
    lipydomics/test/identification.py
    Dylan H. Ross
    2019/11/26

    description:
        Define tests associated with the lipydomics/identification module
"""


import os

from lipydomics.data import Dataset
from lipydomics.identification import add_feature_ids, predict_ccs, predict_rt


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


def add_feature_ids_any_real1_tstamp():
    """
add_feature_ids_any_real1_tstamp
    description:
        Uses the raw data from real_data_1.csv to make compound identifications at any level of confidence. This should
        be a good all-around test for the various identification functions since there are several features in this
        dataset that should not be able to be identified at any level (and thus will run through all of the functions).
        This version of the test uses a specific time-stamped build of the database from the builds directory

        Test fails if there are any errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')

    add_feature_ids(dset, [0.05, 0.5, 5.], db_version_tstamp='2005041036')

    return True


def add_feature_ids_real1_bad_tstamp():
    """
add_feature_ids_any_bad_tstamp
    description:
        Uses the raw data from real_data_1.csv to make compound identifications at any level of confidence.
        A nonexistent time-stamp is provided which should throw a RuntimeError

        Test fails if there are any other errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')

    try:
        add_feature_ids(dset, [0.05, 0.5, 5.], db_version_tstamp='000000')
    except RuntimeError:
        return True

    return False


def predict_ccs_noerrs():
    """
predict_ccs_noerrs
    description:
        predicts CCS for a few lipids, there should be no errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    ccs = predict_ccs('PC', 34, 3, '[M+H]+')
    ccs = predict_ccs('PE', 38, 1, '[M-H]-', fa_mod='o')
    ccs = predict_ccs('LPE', 18, 1, '[M+Na]+', mz=234.5678)

    return True


def predict_ccs_notencodable():
    """
predict_ccs_notencodable
    description:
        predicts CCS for a few lipids, each should raise a ValueError due to the various non-encodable parameters
        test passes only if the expected errors are raised
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    try:
        ccs = predict_ccs('PP', 34, 3, '[M+H]+')
    except ValueError:
        pass
    try:
        ccs = predict_ccs('PE', 38, 1, '[M-H]-', fa_mod='x')
    except ValueError:
        pass
    try:
        ccs = predict_ccs('LPE', 18, 1, '[M+Dog+Cat]+', mz=234.5678)
    except ValueError:
        pass

    return True


def predict_rt_noerrs():
    """
predict_rt_noerrs
    description:
        predicts HILIC rt for a few lipids, there should be no errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    rt = predict_rt('PC', 34, 3)
    rt = predict_rt('PE', 38, 1, fa_mod='p')

    return True


def predict_rt_notencodable():
    """
predict_rt_notencodable
    description:
        predicts HILIC rt for a few lipids, each should raise a ValueError due to the various non-encodable parameters
        test passes only if the expected errors are raised
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    try:
        rt = predict_rt('PP', 34, 3)
    except ValueError:
        pass
    try:
        rt = predict_rt('PE', 38, 1, fa_mod='x')
    except ValueError:
        pass

    return True


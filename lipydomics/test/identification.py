"""
    lipydomics/test/identification.py
    Dylan H. Ross
    2019/11/26

    description:
        Define tests associated with the lipydomics/identification module
"""


import os

from lipydomics.test import run_tests
from lipydomics.data import Dataset
from lipydomics.identification import add_feature_ids, predict_ccs, predict_rt, remove_potential_nonlipids


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


def add_feature_ids_any_ppm_real1():
    """
add_feature_ids_any_ppm_real1
    description:
        Uses the raw data from real_data_1.csv to make compound identifications at any level of confidence using ppm as
        the m/z tolerance. This should
        be a good all-around test for the various identification functions since there are several features in this
        dataset that should not be able to be identified at any level (and thus will run through all of the functions).

        Test fails if there are any errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')

    add_feature_ids(dset, [50., 0.5, 5.], mz_tol_type='ppm')

    return True


def add_feature_ids_custom_real1():
    """
add_feature_ids_custom_real1
    description:
        Uses the raw data from real_data_1.csv to make compound identifications using a custom list of confidence
        levels (pred_mz, meas_rt_ccs, pred_mz_rt_ccs, in reverse order).

        Test fails if there are any errors, or if all three identification levels are not present in the Dataset, or
        if any other identification levels are present in the Dataset
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    add_feature_ids(dset, [0.05, 0.5, 5.], level=['pred_mz_rt_ccs', 'pred_mz_rt', 'pred_mz'])
    found_tm, found_tmr, found_tmrc = False, False, False
    for lvl in dset.feat_id_levels:
        if lvl == 'pred_mz':
            found_tm = True
        elif lvl == 'pred_mz_rt':
            found_tmr = True
        elif lvl == 'pred_mz_rt_ccs':
            found_tmrc = True
        elif lvl != '':
            raise RuntimeError('add_feature_ids_custom_real1: unexpected ID level "{}"'.format(lvl))
    if not found_tm:
        raise RuntimeError('add_feature_ids_custom_real1: did not find ID level "pred_mz" in identifications')
    if not found_tmr:
        raise RuntimeError('add_feature_ids_custom_real1: did not find ID level "pred_mz_rt" in identifications')
    if not found_tmrc:
        raise RuntimeError('add_feature_ids_custom_real1: did not find ID level "pred_mz_rt_ccs" in identifications')
    return True


def add_feature_ids_badcustom_real1():
    """
add_feature_ids_badcustom_real1
    description:
        Uses the raw data from real_data_1.csv to make compound identifications using 2 custom lists of confidence
        levels. Both lists have an invalid ID level in them, one uses "any" while the other has a typo

        Test fails if any errors except for the explicitly caught ones are raised
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    try:
        # the custom list has "any" in it
        add_feature_ids(dset, [0.05, 0.5, 5.], level=['any', 'pred_mz_rt', 'pred_mz'])
    except ValueError:
        pass
    try:
        # the custom list has a typo
        add_feature_ids(dset, [0.05, 0.5, 5.], level=['pred_mz__ccs', 'pred_mz_rt', 'pred_mz'])
    except ValueError:
        pass

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

    add_feature_ids(dset, [0.05, 0.5, 5.], db_version_tstamp='2108240850')

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


def predict_ccs_ignencerr():
    """
predict_ccs_ignencerr
    description:
        predicts CCS for a couple of lipids that are not encodable, but uses the ignore_encoding_errors flag to ignore
        that and predict CCS anyways. The first call should produce a ValueError from LipidMass not being able to
        predict m/z, then the second provides an arbitrary m/z which should not produce any errors

        test passes if the first call raises only the expected error and the second call raises no errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    try:
        ccs = predict_ccs('Rock', 34, 3, '[M+Paper]+', fa_mod='Scissors', ignore_encoding_errors=True)
    except ValueError:
        pass
    ccs = predict_ccs('Rock', 34, 3, '[M+Paper]+', mz=234.5678, fa_mod='Scissors', ignore_encoding_errors=True)

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


def predict_rt_ignencerr():
    """
predict_rt_ignencerr
    description:
        predicts HILIC RT for a couple of lipids that are not encodable, but uses the ignore_encoding_errors flag to
        ignore that and predict RT anyways. There should not be any errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    rt = predict_rt('Elephant', 34, 3, ignore_encoding_errors=True)
    rt = predict_rt('PE', 22, 5, fa_mod='Scissors', ignore_encoding_errors=True)

    return True


def remove_potential_nonlipids_bad_esi_mode():
    """
remove_potential_nonlipids_bad_esi_mode
    description:
        ESI mode of the dataset is not 'pos' or 'neg'
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    try:
        remove_potential_nonlipids(dset)
    except ValueError:
        return True
    return False


def remove_potential_nonlipids_features_not_identified():
    """
remove_potential_nonlipids_features_not_identified
    description:
        tries to remove features without actually running feature identification first
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    try:
        remove_potential_nonlipids(dset)
    except RuntimeError:
        return True
    return False


def remove_potential_nonlipids_features_noerr():
    """
remove_potential_nonlipids_features_noerr
    description:
        performs identifications at 3 levels and each time tries to remove features. There should be no errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    add_feature_ids(dset, [0.05, 0.5, 5.], 'any')
    n_any = remove_potential_nonlipids(dset)
    add_feature_ids(dset, [0.05, 0.5, 5.], 'pred_mz_rt')
    n_pred_mz_rt = remove_potential_nonlipids(dset)
    add_feature_ids(dset, [0.05, 0.5, 5.], 'pred_mz')
    n_pred_mz = remove_potential_nonlipids(dset)
    #print('\nany:', n_any, flush=True)
    #print('pred_mz_rt:', n_pred_mz_rt, flush=True)
    #print('pred_mz:', n_pred_mz, flush=True)
    return True


# references to all of the test functions to be run, and order to run them in
all_tests = [
    add_feature_ids_any_real1,
    add_feature_ids_any_ppm_real1,
    add_feature_ids_custom_real1,
    add_feature_ids_badcustom_real1,
    add_feature_ids_any_real1_tstamp,
    add_feature_ids_real1_bad_tstamp,
    predict_ccs_noerrs,
    predict_ccs_notencodable,
    predict_ccs_ignencerr,
    predict_rt_noerrs,
    predict_rt_notencodable,
    predict_rt_ignencerr,
    remove_potential_nonlipids_bad_esi_mode,
    remove_potential_nonlipids_features_not_identified,
    remove_potential_nonlipids_features_noerr
]
if __name__ == '__main__':
    run_tests(all_tests)

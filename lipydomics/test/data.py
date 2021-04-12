"""
    lipydomics/test/data.py
    Dylan H. Ross
    2019/09/25

    description:
        Define tests associated with the lipydomics/data.py module
"""


import os
import numpy as np

from lipydomics.test import run_tests
from lipydomics.data import Dataset
from lipydomics.stats import add_pca3, add_plsra
from lipydomics.identification import add_feature_ids


def dataset_init_mock1():
    """
dataset_init_mock1
    description:
        Tests initialization of the Dataset class using mock_data_1.csv:
            1 header line to skip
            5 features
            6 samples
        --> labels should have shape: (5, 3)
        --> intensities should have shape: (5, 6)
        Test fails if these shapes are incorrect
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    # initialize Dataset
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    return dset.labels.shape == (5, 3) and dset.intensities.shape == (5, 6)


def dataset_normalize_mock1():
    """
dataset_normalize_mock1
    description:
        TODO
        Test fails if
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    weights = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.normalize(weights)
    return True


def dataset_getgroup_mock1():
    """
dataset_getgroup_mock1
    description:
        Assigns groups "A" and "B" to the dataset with some indices, then gets the group data by group name
        Test fails only if there are errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({"A": [0, 2, 4], "B": [1, 3, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    A, B = dset.get_data_bygroup(["A", "B"], normed=True)
    return True


def dataset_assign_groups_using_replicates_real1():
    """
dataset_assign_groups_using_replicates_real1
    description:
        Loads the real data example from .csv then uses Dataset.assign_groups_and_replicates to assign groups and checks
        that the indices are correct.

        Test fails if Dataset.group_indices is not set as expected.
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups_with_replicates(['Par', 'Dap2', 'Dal2', 'Van4', 'Van8'], 4)
    group_indices = {
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7],
        'Dal2': [8, 9, 10, 11],
        'Van4': [12, 13, 14, 15],
        'Van8': [16, 17, 18, 19]
    }
    for group in dset.group_indices.keys():
        if dset.group_indices[group] != group_indices[group]:
            # indices are not the same
            return False
    # no group index mismatches were found
    return True


def dataset_save_load_bin_mock1():
    """
dataset_save_load_bin_mock1
    description:
        Loads the mock data from .csv, assigns groups, normalizes the data then saves the instance to serialized binary
        format. The binary file is then reloaded and presence of the groups and normalized data is verified.

        Test fails if there are any errors, if the .pickle file does not get produced, or if the reloaded Dataset does
        not have the correct attributes.
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({"A": [0, 2, 4], "B": [1, 3, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    bin_path = os.path.join(os.path.dirname(__file__), 'dataset_mock1.pickle')
    dset.save_bin(bin_path)
    if not os.path.isfile(bin_path):
        m = 'dataset_save_load_bin_mock1: binary file {} not found'
        raise RuntimeError(m.format(bin_path))
    # re-load the dataset
    dset_reload = Dataset.load_bin(bin_path)
    if dset.group_indices is None:
        raise RuntimeError('dataset_save_load_bin_mock1: dset_reload.group_indices not set')
    if dset.normed_intensities is None:
        raise RuntimeError('dataset_save_load_bin_mock1: dset_reload.normed_intensities not set')
    # remove the binary file
    os.remove(bin_path)
    return True


def dataset_export_feature_data_real1():
    """
dataset_export_feature_data_real1
    description:
        Loads a dataset (real_data_1.csv) then exports a selection of features, first with a Dataset that has only
        raw data, then with a Dataset that has raw and normalized data.

        Test fails if there are any errors, or if either of the output .csv files do not get produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    incsv = os.path.join(os.path.dirname(__file__), 'real1_features1.csv')
    outcsv1 = os.path.join(os.path.dirname(__file__), 'real1_export_raw1.csv')
    outcsv2 = os.path.join(os.path.dirname(__file__), 'real1_export_norm1.csv')
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.select_feature_data(incsv, outcsv1)
    # make sure the output file gets created
    if not os.path.isfile(outcsv1):
        m = 'dataset_export_feature_data_real1: output {} not found'
        raise RuntimeError(m.format(outcsv1))
    # normalize the data and try again
    dset.normalize(np.array([0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75,
                             0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75]))
    dset.select_feature_data(incsv, outcsv2)
    if not os.path.isfile(outcsv2):
        m = 'dataset_export_feature_data_real1: output {} not found'
        raise RuntimeError(m.format(outcsv2))
    # remove the output files
    os.remove(outcsv1)
    os.remove(outcsv2)
    return True


def dataset_export_xlsx_real1():
    """
dataset_export_xlsx_real1
    description:
        Loads a dataset (real_data_1.csv) then exports the dataset as an Excel spreadsheet

        Test fails if there are any errors, or if the spreadsheet is not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    out = os.path.join(os.path.dirname(__file__), 'real1_exported.xlsx')
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    # normalize the data
    dset.normalize(np.array([0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75,
                             0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75]))
    dset.export_xlsx(out)
    # make sure the output file gets created
    if not os.path.isfile(out):
        m = 'dataset_export_xlsx_real1: output {} not found'
        raise RuntimeError(m.format(out))
    # remove the exported file
    os.remove(out)
    return True


def dataset_export_analyzed_xlsx_real1():
    """
dataset_export_analyzed_xlsx_real1
    description:
        Loads a dataset (real_data_1.csv) then performs some analyses and exports the dataset as an Excel spreadsheet

        Test fails if there are any errors, or if the spreadsheet is not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    out = os.path.join(os.path.dirname(__file__), 'real1_exported.xlsx')
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    # normalize the data
    dset.normalize(np.array([0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75,
                             0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75]))
    # compute a PCA then identify lipids without using RT
    dset.assign_groups_with_replicates(['Par', 'Dap2', 'Dal2', 'Van4', 'Van8'], 4)
    add_pca3(dset, ['Par', 'Dap2', 'Dal2', 'Van4', 'Van8'], normed=True)
    add_feature_ids(dset, [0.025, 1.0, 3.0], level='any', use_rt=False)
    dset.export_xlsx(out)
    # make sure the output file gets created
    if not os.path.isfile(out):
        m = 'dataset_export_xlsx_real1: output {} not found'
        raise RuntimeError(m.format(out))
    # remove the exported file
    os.remove(out)
    return True


def dataset_drop_features_goodparams_real1():
    """
dataset_drop_features_goodparams_real1
    description:
        Loads a dataset (real_data_1.csv) then calls the drop_features(...) with a few parameter combinations

        Test fails if any errors occur or if any of the resulting arrays are not of the expected shape
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    np.random.seed(69)
    stat1 = np.random.random((773,))
    stat3 = np.random.random((3, 773))
    stat3t = np.random.random((773, 3))
    try_params = [
        ('mintensity', {'normed': False, 'lower_bound': 1000}),
        ('meantensity', {'normed': False, 'lower_bound': 50}),
        ('meantensity', {'normed': False, 'lower_bound': 50, 'upper_bound': 5000}),
        ('mock_stat_1column', {'lower_bound': 0.1}),
        ('mock_stat_1column', {'lower_bound': 0.1, 'upper_bound': 0.9}),
        ('mock_stat_3column', {'lower_bound': 0.1, 'upper_bound': 0.9, 'axis': 0}),
        ('mock_stat_3column', {'lower_bound': 0.1, 'upper_bound': 0.9, 'axis': 2}),
        ('mock_stat_3column_T', {'lower_bound': 0.1, 'upper_bound': 0.9, 'axis': 0}),
        ('mock_stat_3column_T', {'lower_bound': 0.1, 'upper_bound': 0.9, 'axis': 2}),
    ]
    for c, kw in try_params:
        # load the data
        dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
        # add in some mock statistics
        dset.stats['mock_stat_1column'] = stat1
        dset.stats['mock_stat_3column'] = stat3
        dset.stats['mock_stat_3column_T'] = stat3t
        # drop features
        dset.drop_features(c, **kw)
        # validate that the drop worked by checking the resulting shapes
        assert dset.n_features < 773
        assert dset.intensities.shape == (dset.n_features, dset.n_samples)
        assert dset.labels.shape == (dset.n_features, 3)
        assert dset.stats['mock_stat_1column'].shape == (dset.n_features,)
        assert dset.stats['mock_stat_3column'].shape == (3, dset.n_features)
        assert dset.stats['mock_stat_3column_T'].shape == (dset.n_features, 3)
    return True


def dataset_drop_features_badparams_real1():
    """
dataset_drop_features_badparams_real1
    description:
        Loads a dataset (real_data_1.csv) then calls the drop_features(...) method with a bunch of bad parameters meant
        to trigger ValueErrors

        Test fails if any uncaught errors occur
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    # add in some mock statistics
    dset.stats['mock_stat_1column'] = np.random.random((dset.n_features,))
    dset.stats['mock_stat_3column'] = np.random.random((3, dset.n_features))
    try_params = [
        ('kitten', {}),
        ('meantensity', {}),
        ('mintensity', {'normed': False}),
        ('mintensity', {'normed': False, 'upper_bound': 100}),
        ('meantensity', {'normed': True, 'lower_bound': 100}),
        ('mock_stat_1column', {}),
        ('mock_stat_3column', {'lower_bound': 100}),
    ]
    for c, kw in try_params:
        try:
            dset.drop_features(c, **kw)
        except ValueError as ve:
            #print(ve)
            pass
    return True


def dataset_drop_features_real1():
    """
dataset_drop_features_real1
    description:
        Loads a dataset (real_data_1.csv), works up the data by identifying lipids and computing statistics, then calls
        the drop_features(...) method (using meantensity with two bounds)

        Test fails if any errors occur or if any of the components of the Dataset are not the expected shape after the
        transformation
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    dset.assign_groups_with_replicates(['A', 'B', 'C', 'D', 'E'], 4)
    add_pca3(dset, ['A', 'B', 'C', 'D', 'E'])
    y = []
    for i in range(5):
        for j in range(4):
            y.append(float(i))
    y = np.array(y)
    add_plsra(dset, ['A', 'B', 'C', 'D', 'E'], y)
    add_feature_ids(dset, [0.05, 0.2, 0.1], level='pred_mz_rt')
    #print(dset)
    dset.drop_features('meantensity', normed=False, lower_bound=50, upper_bound=1000)
    #print(dset)
    # validate that the dropping worked by checking some shapes
    assert dset.n_features < 773
    assert dset.intensities.shape == (dset.n_features, dset.n_samples)
    assert dset.labels.shape == (dset.n_features, 3)
    assert dset.stats['PCA3_A-B-C-D-E_loadings_raw'].shape == (3, dset.n_features)
    assert dset.stats['PCA3_A-B-C-D-E_projections_raw'].shape == (dset.n_samples, 3)
    assert dset.ext_var.shape == (dset.n_samples,)
    assert dset.stats['PLS-RA_A-B-C-D-E_loadings_raw'].shape == (dset.n_features, 2)
    assert dset.stats['PLS-RA_A-B-C-D-E_projections_raw'].shape == (dset.n_samples, 2)
    assert len(dset.feat_ids) == dset.n_features
    assert len(dset.feat_id_levels) == dset.n_features
    assert len(dset.feat_id_scores) == dset.n_features
    return True


# references to al of the test functions to be run, and order to run them in
all_tests = [
    dataset_init_mock1,
    dataset_normalize_mock1,
    dataset_getgroup_mock1,
    dataset_assign_groups_using_replicates_real1,
    dataset_save_load_bin_mock1,
    dataset_export_feature_data_real1,
    dataset_export_xlsx_real1,
    dataset_export_analyzed_xlsx_real1,
    dataset_drop_features_goodparams_real1,
    dataset_drop_features_badparams_real1,
    dataset_drop_features_real1
]
if __name__ == '__main__':
    run_tests(all_tests)

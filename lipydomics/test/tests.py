"""
    lipydomics/test/tests.py
    Dylan H. Ross
    2019/02/02

    description:
        TODO
"""


import traceback
import os
import numpy as np

from lipydomics.data import Dataset
from lipydomics.stats import add_anova_p


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


def stats_addanovap_mock1():
    """
stats_addanovap_mock1
    description:
        Uses the normalized data from mock_data_1.csv to compute ANOVA p-values for groups "A" and "B"
        Test fails if the shape of the resulting "test_anova" statistic array in the Dataset object does not have 
        the proper shape: (5,)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({"A": [0, 2, 4], "B": [1, 3, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    add_anova_p(dset, ["A", "B"], normed=True, label="test_anova")
    return dset.stats["test_anova"].shape == (5,)


def stats_addpca3loadings_mock1():
    """
stats_addpca3loadings_mock1
    description:
       
       Test fails if
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    return False


def run_all_tests():
    """
run_all_tests
    description:
        runs all tests sequentially, if there are any failures the test function docstring is printed
"""
    # references to al of the test functions to be run
    all_tests = [
        dataset_init_mock1,
        dataset_normalize_mock1,
        dataset_getgroup_mock1,
        stats_addanovap_mock1,
        stats_addpca3loadings_mock1
    ]
    # run the tests
    failed = False
    print("running all tests ... ", end="")
    for test in all_tests:
        try:
            passed = test()
        except Exception as e:
            print('\n', test.__doc__, 'TEST FAILED WITH EXEPTION!\n', e)
            print(traceback.format_exc())
            failed = True
            break
        if not passed:
            print('\n', test.__doc__, 'TEST FAILED!\n')
            failed = True
            break
    if not failed:
        print("passed")


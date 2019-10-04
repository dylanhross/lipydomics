"""
    lipydomics/test/data.py
    Dylan H. Ross
    2019/09/25

    description:
        Define tests associated with the lipydomics/data.py module
"""


import os
import numpy as np

from lipydomics.data import Dataset


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


def dataset_save_load_bin_mock1():
    """
dataset_save_load_bin_mock1
    description:
        Loads the mock data from .csv, assigns groups, normalizes the data then saves the instance to serialized binary
        format. The binary file is then reloaded and presence of the groups and normalized data is verified.

        Test fails if there are any errors, if the .pickle file does not get produced, or if the reloaded Dataset does
        not have the coorrect attributes. 
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


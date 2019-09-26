"""
    lipydomics/test/stats.py
    Dylan H. Ross
    2019/09/25

    description:
        Define tests associated with the lipydomics/stats.py module
"""


import os
import numpy as np

from lipydomics.data import Dataset
from lipydomics.stats import add_anova_p, add_pca3, add_plsda


def addanovap_mock1():
    """
addanovap_mock1
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
    add_anova_p(dset, ["A", "B"], normed=True)
    return dset.stats["ANOVA_A-B_normed"].shape == (5,)


def addpca3_mock1():
    """
addpca3_mock1
    description:
        Uses the raw data from mock_data_1.csv to compute a 3 component PCA.

        Test fails if there are any errors or if the shapes of the following stats entries are incorrect:
            dset.stats['PCA3_loadings_raw'].shape = (3, 5)
            dset.stats['PCA3_projections_raw'].shape = (6, 3)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 2, 4], 'B': [1, 3, 5]})
    add_pca3(dset)
    if dset.stats['PCA3_loadings_raw'].shape != (3, 5):
        m = 'addpca3_mock1: PCA3_loadings_raw should have shape (3, 5), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PCA3_loadings_raw'].shape))
    if dset.stats['PCA3_projections_raw'].shape != (6, 3):
        m = 'addpca3_mock1: PCA3_projections_raw should have shape (6, 3), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PCA3_projections_raw'].shape))

    return True


def addplsda_mock1():
    """
addplsda_mock1
    description:
        Uses the raw data from mock_data_1.csv to perform PLS-DA.

        Test fails if there are any errors or if the shapes of the following stats entries are incorrect:
            dset.stats['PLS-DA_A-B_loadings_raw'].shape = (5, 2)
            dset.stats['PLS-DA_A-B_projections_raw'].shape = (6, 2)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 2, 4], 'B': [1, 3, 5]})
    add_plsda(dset, ['A', 'B'])
    if dset.stats['PLS-DA_A-B_loadings_raw'].shape != (5, 2):
        m = 'addplsda_mock1: PCA3_loadings_raw should have shape (5, 2), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PLS-DA_A-B_loadings_raw'].shape))
    if dset.stats['PLS-DA_A-B_projections_raw'].shape != (6, 2):
        m = 'addplsda_mock1: PCA3_projections_raw should have shape (6, 2), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PLS-DA_A-B_projections_raw'].shape))

    return True



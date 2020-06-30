"""
    lipydomics/test/stats.py
    Dylan H. Ross
    2019/09/25

    description:
        Define tests associated with the lipydomics/stats.py module
"""


import os
import numpy as np

from lipydomics.test import run_tests
from lipydomics.data import Dataset
from lipydomics.stats import (
    add_anova_p, add_pca3, add_plsda, add_2group_corr, add_plsra, add_log2fc, add_2group_pvalue
)


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


def addanovap_real1():
    """
addanovap_real1
    description:
        Uses the raw data from real_data_1.csv to compute ANOVA p-values for all 5 of the defined groups
        Test fails if the shape of the resulting "test_anova" statistic array in the Dataset object does not have 
        the proper shape: (773,)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3], 
        'Dap2': [4, 5, 6, 7], 
        'Dal2': [8, 9, 10, 11], 
        'Van4': [12, 13, 14, 15], 
        'Van8': [16, 17, 18, 19]
    })
    group_names = ['Par', 'Dap2', 'Dal2', 'Van4', 'Van8']
    add_anova_p(dset, group_names)
    return dset.stats['ANOVA_{}_raw'.format('-'.join(group_names))].shape == (773,)


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
    add_pca3(dset, ['A', 'B'])
    if dset.stats['PCA3_A-B_loadings_raw'].shape != (3, 5):
        m = 'addpca3_mock1: PCA3_A-B_loadings_raw should have shape (3, 5), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PCA3_A-B_loadings_raw'].shape))
    if dset.stats['PCA3_A-B_projections_raw'].shape != (6, 3):
        m = 'addpca3_mock1: PCA3_A-B_projections_raw should have shape (6, 3), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PCA3_A-B_projections_raw'].shape))

    return True


def addpca3_real1():
    """
addpca3_real1
    description:
        Uses the raw data from real_data_1.csv to compute a 3 component PCA.

        Test fails if there are any errors or if the shapes of the following stats entries are incorrect:
            dset.stats['PCA3_..._loadings_raw'].shape = (3, 773)
            dset.stats['PCA3_..._projections_raw'].shape = (20, 3)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3], 
        'Dap2': [4, 5, 6, 7], 
        'Dal2': [8, 9, 10, 11], 
        'Van4': [12, 13, 14, 15], 
        'Van8': [16, 17, 18, 19]
    })
    group_names = ['Par', 'Dap2', 'Dal2', 'Van4', 'Van8']
    add_pca3(dset, group_names)
    if dset.stats['PCA3_{}_loadings_raw'.format('-'.join(group_names))].shape != (3, 773):
        m = 'addpca3_real1: PCA3_{}_loadings_raw should have shape (3, 773), has shape: {}'
        raise RuntimeError(m.format('-'.join(group_names), dset.stats['PCA3_{}_loadings_raw'.format('-'.join(group_names))].shape))
    if dset.stats['PCA3_{}_projections_raw'.format('-'.join(group_names))].shape != (20, 3):
        m = 'addpca3_real1: PCA3_{}_projections_raw should have shape (20, 3), has shape: {}'
        raise RuntimeError(m.format('-'.join(group_names), dset.stats['PCA3_{}_projections_raw'.format('-'.join(group_names))].shape))

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
        m = 'addplsda_mock1: PLS-DA_A-B_loadings_raw should have shape (5, 2), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PLS-DA_A-B_loadings_raw'].shape))
    if dset.stats['PLS-DA_A-B_projections_raw'].shape != (6, 2):
        m = 'addplsda_mock1: PLS-DA_A-B_projections_raw should have shape (6, 2), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PLS-DA_A-B_projections_raw'].shape))

    return True


def addplsda_3groups_mock1():
    """
addplsda_3groups_mock1
    description:
        Uses the raw data from mock_data_1.csv to perform PLS-DA trying to use 3 groups

        Test fails if there are any unexpected errors or if the expected ValueError does not occur
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 2, 4], 'B': [1, 3, 5]})
    try:
        add_plsda(dset, ['A', 'B', 'C'])
    except ValueError:
        return True
    
    return False


def addplsda_real1():
    """
addplsda_real1
    description:
        Uses the raw data from real_data_1.csv to perform PLS-DA.

        Test fails if there are any errors or if the shapes of the following stats entries are incorrect:
            dset.stats['PLS-DA_..._loadings_raw'].shape = (773, 2)
            dset.stats['PLS-DA_..._projections_raw'].shape = (8, 2)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3], 
        'Dap2': [4, 5, 6, 7], 
        'Dal2': [8, 9, 10, 11], 
        'Van4': [12, 13, 14, 15], 
        'Van8': [16, 17, 18, 19]
    })

    # pairs of groups to compute PLS-DA on
    pairs = [
        ['Par', 'Dap2'],
        ['Par', 'Dal2'],
        ['Par', 'Van4'],
        ['Par', 'Van8']
    ]

    for pair in pairs:
        add_plsda(dset, pair)
        if dset.stats['PLS-DA_{}_loadings_raw'.format('-'.join(pair))].shape != (773, 2):
            m = 'addplsda_real1: PLS-DA_loadings_raw should have shape (773, 2), has shape: {}'
            raise RuntimeError(m.format(dset.stats['PLS-DA_{}_loadings_raw'.format('-'.join(pair))].shape))
        if dset.stats['PLS-DA_{}_projections_raw'.format('-'.join(pair))].shape != (8, 2):
            m = 'addplsda_real1: PLS-DA_projections_raw should have shape (8, 2), has shape: {}'
            raise RuntimeError(m.format(dset.stats['PLS-DA_{}_projections_raw'.format('-'.join(pair))].shape))

    return True


def add2groupcorr_mock1():
    """
add2groupcorr_mock1
    description:
        Uses the normalized data from mock_data_1.csv to compute Pearson correlation for groups "A" and "B"
        Test fails if the shape of the resulting statistic array in Dataset.stats does not have 
        the proper shape: (5,)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({"A": [0, 2, 4], "B": [1, 3, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    add_2group_corr(dset, ["A", "B"], normed=True)
    return dset.stats['2-group-corr_A-B_normed'].shape == (5,)


def add2groupcorr_3groups_mock1():
    """
add2groupcorr_3groups_mock1
    description:
        Uses the raw data from mock_data_1.csv to try to compute Pearson correlation trying to use 3 groups

        Test fails if there are any unexpected errors or if the expected ValueError does not occur
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({"A": [0, 2, 4], "B": [1, 3, 5], "C": [1, 2, 3]})
    try:
        add_2group_corr(dset, ['A', 'B', 'C'])
    except ValueError:
        return True
    
    return False


def add2groupcorr_real1():
    """
add2groupcorr_real1
    description:
        Uses the raw data from real_data_1.csv to compute Pearson correlation for groups "A" and "B"
        Test fails if the shape of the resulting statistic array in Dataset.stats does not have 
        the proper shape: (5,)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3], 
        'Dap2': [4, 5, 6, 7], 
        'Dal2': [8, 9, 10, 11], 
        'Van4': [12, 13, 14, 15], 
        'Van8': [16, 17, 18, 19]
    })

    # pairs of groups to compute Pearson correlation on
    pairs = [
        ['Par', 'Dap2'],
        ['Par', 'Dal2'],
        ['Par', 'Van4'],
        ['Par', 'Van8']
    ]

    for pair in pairs:
        add_2group_corr(dset, pair)
        if dset.stats['2-group-corr_{}_raw'.format('-'.join(pair))].shape != (773,):
            m = 'add2groupcorr_real1: PLS-DA_loadings_raw should have shape (773, 2), has shape: {}'
            raise RuntimeError(m.format(dset.stats['2-group-corr_{}_raw'.format('-'.join(pair))].shape))

    return True


def addplsra_real1():
    """
addplsra_real1
    description:
        Uses the raw data from real_data_1.csv to perform PLS-RA with all groups.

        Test fails if there are any errors or if the shapes of the following stats entries are incorrect:
            dset.stats['PLS-DA_..._loadings_raw'].shape = (773, 2)
            dset.stats['PLS-DA_..._projections_raw'].shape = (20, 2)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7],
        'Dal2': [8, 9, 10, 11],
        'Van4': [12, 13, 14, 15],
        'Van8': [16, 17, 18, 19]
    })

    with open(os.path.join(os.path.dirname(__file__), 'external_real1.txt'), 'r') as f:
        y = np.array([float(_.strip()) for _ in f.readlines()])

    groups = ['Par', 'Dap2', 'Dal2', 'Van4', 'Van8']
    add_plsra(dset, groups, y)
    if dset.stats['PLS-RA_{}_loadings_raw'.format('-'.join(groups))].shape != (773, 2):
        m = 'addplsra_real1: PLS-RA_loadings_raw should have shape (773, 2), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PLS-RA_{}_loadings_raw'.format('-'.join(groups))].shape))
    if dset.stats['PLS-RA_{}_projections_raw'.format('-'.join(groups))].shape != (20, 2):
        m = 'addplsra_real1: PLS-RA_projections_raw should have shape (20, 2), has shape: {}'
        raise RuntimeError(m.format(dset.stats['PLS-RA_{}_projections_raw'.format('-'.join(groups))].shape))

    return True


def addlog2fc_real1():
    """
addlog2fc_real1
    description:
        Uses the raw data from real_data_1.csv to compute log2(fold-change) on a bunch of pairs of groups

        Test fails if there are any errors or if the shape of any of the following stats entries are incorrect:
            dset.stats['log2fc_..._raw'].shape = (773,)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7],
        'Dal2': [8, 9, 10, 11],
        'Van4': [12, 13, 14, 15],
        'Van8': [16, 17, 18, 19]
    })

    # pairs of groups to compute log2fc on
    pairs = [
        ['Par', 'Dap2'],
        ['Par', 'Dal2'],
        ['Par', 'Van4'],
        ['Par', 'Van8']
    ]

    for pair in pairs:
        add_log2fc(dset, pair)
        if dset.stats['LOG2FC_{}_raw'.format('-'.join(pair))].shape != (773,):
            m = 'addlog2fc_real1: "log2fc_..._raw" should have shape (773,), has shape: {}'
            raise RuntimeError(m.format(dset.stats['LOG2FC_{}_raw'.format('-'.join(pair))].shape))

    return True


def add2grouppvalue_real1():
    """
add2grouppvalue_real1
    description:
        Uses the raw data from real_data_1.csv to compute p-values for 2 group comparisons on a bunch of pairs of
        groups and using all 3 stats tests

        Test fails if there are any errors or if the shape of any of the following stats entries are incorrect:
            dset.stats['{stats_test}_..._raw'].shape = (773,)
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7],
        'Dal2': [8, 9, 10, 11],
        'Van4': [12, 13, 14, 15],
        'Van8': [16, 17, 18, 19]
    })

    # pairs of groups to compute log2fc on
    pairs = [
        ['Par', 'Dap2'],
        ['Par', 'Dal2'],
        ['Par', 'Van4'],
        ['Par', 'Van8']
    ]

    for pair in pairs:
        for stats_test in ['students', 'welchs', 'mann-whitney']:
            #print('testing {} with {}'.format(pair, stats_test))
            add_2group_pvalue(dset, pair, stats_test)
            stest_abbrev = {'students': 'studentsP', 'welchs': 'welchsP', 'mann-whitney': 'mannwhitP'}[stats_test]
            if dset.stats['{}_{}_raw'.format(stest_abbrev, '-'.join(pair))].shape != (773,):
                m = 'add2grouppvalue_real1: "{}_..._raw" should have shape (773,), has shape: {}'
                raise RuntimeError(m.format(dset.stats['LOG2FC_{}_raw'.format(stest_abbrev, '-'.join(pair))].shape))

    # diagnostic printing stuff
    """
    print(dset)
    for s, w, m in zip(dset.stats["studentsP_Par-Dap2_raw"] <= 0.05,
                       dset.stats["welchsP_Par-Dap2_raw"] <= 0.05,
                       dset.stats["mannwhitP_Par-Dap2_raw"] <= 0.05):
        if s and w and m:
            print(True)
        elif not s and not w and not m:
            print(False)
        else:
            print(s, w, m)
    """

    return True


# references to al of the test functions to be run, and order to run them in
all_tests = [
    addanovap_mock1,
    addanovap_real1,
    addpca3_mock1,
    addpca3_real1,
    addplsda_mock1,
    addplsda_3groups_mock1,
    addplsda_real1,
    add2groupcorr_mock1,
    add2groupcorr_3groups_mock1,
    add2groupcorr_real1,
    addplsra_real1,
    addlog2fc_real1,
    add2grouppvalue_real1
]
if __name__ == '__main__':
    run_tests(all_tests)

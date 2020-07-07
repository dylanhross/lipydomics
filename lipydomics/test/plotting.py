"""
    lipydomics/test/plotting.py
    Dylan H. Ross
    2019/09/25

    description:
        Define tests associated with the lipydomics/plotting.py module
"""


import os
import numpy as np
from csv import reader

from lipydomics.test import run_tests
from lipydomics.data import Dataset
from lipydomics.stats import add_pca3, add_plsda, add_2group_corr, add_plsra, add_log2fc, add_2group_pvalue
from lipydomics.plotting import (
    barplot_feature_bygroup, batch_barplot_feature_bygroup, scatter_pca3_projections_bygroup,
    scatter_plsda_projections_bygroup, splot_plsda_pcorr_bygroup, scatter_plsra_projections_bygroup,
    heatmap_lipid_class_log2fc, volcano_2group
)
from lipydomics.identification import add_feature_ids


def barplot_feature_bygroup_mock1():
    """
barplot_feature_bygroup_mock1
    description:
        Using the normalized data from mock_data_1.csv, extract and plot average values for each feature. 
        Image files are saved in the lipydomics/test directory, then deleted. 

        Test fails if there are any errors, or if the image files are not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 2, 4], 'B': [1, 3, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    
    png_path = os.path.dirname(__file__)
    features = [
        (123.4567, 1.23, 200.46),
        (234.5678, 2.34, 244.68),
        (345.6789, 3.45, 266.80),
        (456.7890, 4.56, 288.02),
        (567.8901, 5.67, 300.24)
    ]

    # process each feature, then check for the image file 
    group_names = ['A', 'B']
    for feat in features:
        img_dir = os.path.dirname(__file__)
        img_path = os.path.join(img_dir, 'bar_{:.4f}-{:.2f}-{:.1f}_{}_normed.png'.format(*feat, '-'.join(group_names)))
        barplot_feature_bygroup(dset, group_names, img_dir, feat, normed=True)
        if not os.path.isfile(img_path):
            m = 'barplot_feature_bygroup_mock1: image file {} not found'
            raise RuntimeError(m.format(img_path))
        else:
            # delete the image file 
            os.remove(img_path)

    return True


def batch_barplot_feature_bygroup_real1():
    """
batch_barplot_feature_bygroup_mock1
    description:
        Using the normalized data from real_data_1.csv, extract and plot average values for each feature defined in
        real1_features1.csv
        Image files are saved in the lipydomics/test directory, then deleted.

        Test fails if there are any errors, or if the image files are not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'))
    dset.assign_groups({'A': [0, 2, 4], 'B': [1, 3, 5]})
    dset.normalize(np.array([0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75,
                             0.75, 0.8, 0.825, 0.85, 0.95, 0.95, 0.85, 0.825, 0.8, 0.75]))
    group_names = ['Par', 'Dap2', 'Dal2', 'Van4', 'Van8']
    dset.assign_groups_with_replicates(group_names, 4)
    img_dir = os.path.dirname(__file__)
    in_csv = os.path.join(os.path.dirname(__file__), 'real1_features1.csv')
    batch_barplot_feature_bygroup(dset, group_names, img_dir, in_csv, normed=True)
    # go through and remove all of the barplots
    with open(in_csv, 'r') as inf:
        next(inf)
        rdr = reader(inf)
        for mz, rt, ccs in rdr:
            mz, rt, ccs = float(mz), float(rt), float(ccs)
            img_path = os.path.join(img_dir,
                                    'bar_{:.4f}-{:.2f}-{:.1f}_{}_normed.png'.format(mz, rt, ccs, '-'.join(group_names)))
            if not os.path.isfile(img_path):
                m = 'batch_barplot_feature_bygroup_real1: image file {} not found'
                raise RuntimeError(m.format(img_path))
            else:
                # delete the image file
                os.remove(img_path)
    return True


def scatter_pca3_projections_bygroup_mock1():
    """
scatter_pca3_projections_bygroup_mock1
    description:
        Using the normalized data from mock_data_1.csv, perform a 3-component PCA and plot the projections 
        Image file is saved in the lipydomics/test directory, then deleted. 

        Test fails if there are any errors, or if the image file is not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 1], 'B': [2, 3], 'C': [4, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    group_names = ['A', 'B', 'C']
    add_pca3(dset, group_names, normed=True)

    img_dir = os.path.dirname(__file__)
    img_path = os.path.join(img_dir, 'PCA3_{}_projections_normed.png'.format('-'.join(group_names)))

    scatter_pca3_projections_bygroup(dset, group_names, img_dir, normed=True)
    if not os.path.isfile(img_path):
        m = 'scatter_pca3_projections_bygroup_mock1: image file {} not found'
        raise RuntimeError(m.format(img_path))
    else:
        # delete the image file 
        os.remove(img_path)

    return True


def scatter_plsda_projections_bygroup_mock1():
    """
scatter_plsda_projections_bygroup_mock1
    description:
        Using the normalized data from mock_data_1.csv, perform PLS-DA and plot the projections 
        Image file is saved in the lipydomics/test directory, then deleted. 

        Test fails if there are any errors, or if the image file is not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 1, 2], 'B': [3, 4, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    group_names = ['A', 'B']
    add_plsda(dset, group_names, normed=True)

    img_dir = os.path.dirname(__file__)
    img_path = os.path.join(img_dir, 'PLS-DA_{}_projections_normed.png'.format('-'.join(group_names)))

    scatter_plsda_projections_bygroup(dset, group_names, img_dir, normed=True)
    if not os.path.isfile(img_path):
        m = 'scatter_plsda_projections_bygroup_mock1: image file {} not found'
        raise RuntimeError(m.format(img_path))
    else:
        # delete the image file 
        os.remove(img_path)

    return True


def splot_plsda_pcorr_bygroup_mock1():
    """
splot_plsda_pcorr_bygroup_mock1
    description:
        Using the normalized data from mock_data_1.csv, perform PLS-DA and correlation analyses then generate an S-plot 
        Image file is saved in the lipydomics/test directory, then deleted. 

        Test fails if there are any errors, or if the image file is not produced
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'mock_data_1.csv'))
    dset.assign_groups({'A': [0, 1, 2], 'B': [3, 4, 5]})
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]))
    group_names = ['A', 'B']
    add_plsda(dset, group_names, normed=True)
    add_2group_corr(dset, group_names, normed=True)

    img_dir = os.path.dirname(__file__)
    img_path = os.path.join(img_dir, 'S-Plot_{}_normed.png'.format('-'.join(group_names)))

    splot_plsda_pcorr_bygroup(dset, group_names, img_dir, normed=True)
    if not os.path.isfile(img_path):
        m = 'scatter_plsda_projections_bygroup_mock1: image file {} not found'
        raise RuntimeError(m.format(img_path))
    else:
        # delete the image file 
        os.remove(img_path)
        
    return True


def scatter_plsra_projections_bygroup_real1():
    """
scatter_plsra_projections_bygroup_real1
    description:
        Using the normalized data from real_data_1.csv, perform PLS-RA and plot the projections
        Image file is saved in the lipydomics/test directory, then deleted.

        Test fails if there are any errors, or if the image file is not produced
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
    dset.normalize(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.1, 0.2, 0.3, 0.4,
                             0.5, 0.6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.1, 0.2]))
    with open(os.path.join(os.path.dirname(__file__), 'external_real1.txt'), 'r') as f:
        y = np.array([float(_.strip()) for _ in f.readlines()])
    group_names = ['Par', 'Dap2', 'Dal2', 'Van4', 'Van8']
    add_plsra(dset, group_names, y, normed=True)

    img_dir = os.path.dirname(__file__)
    img_path = os.path.join(img_dir, 'PLS-RA_{}_projections_normed.png'.format('-'.join(group_names)))

    scatter_plsra_projections_bygroup(dset, group_names, img_dir, normed=True)
    if not os.path.isfile(img_path):
        m = 'scatter_plsra_projections_bygroup_real1: image file {} not found'
        raise RuntimeError(m.format(img_path))
    else:
        # delete the image file
        os.remove(img_path)

    return True


def heatmap_lipid_class_log2fc_real1():
    """
fetch_lipid_class_log2fa_real1
    description:
        Tests the function that makes a heatmap of lipid class fold-change data using identifications made on the
        real_data_1.csv test dataset

        Test fails if there are any errors, if the data is not  found, or the image file is not created
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    # setup
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    dset.assign_groups({
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7]
    })
    add_log2fc(dset, ['Par', 'Dap2'])
    add_feature_ids(dset, [0.05, 0.5, 3.0])

    # plot the PGs
    img_dir = os.path.dirname(__file__)
    fig_path = os.path.join(img_dir, 'PG_Par-Dap2_log2fc_raw.png')
    if not heatmap_lipid_class_log2fc('PG', dset, ['Par', 'Dap2'], img_dir):
        m = 'heatmap_lipid_class_log2fc_real1: no PG data found in identifications'
        raise RuntimeError(m)

    if not os.path.isfile(fig_path):
        m = 'heatmap_lipid_class_log2fc_real1: image file {} not found'
        raise RuntimeError(m.format(fig_path))
    else:
        # delete the image file
        os.remove(fig_path)

    return True


def volcano_2group_real1():
    """
volcano_2group_real1
    description:
        Tests the function that makes a volcano plot

        Test fails if there are any errors or the image file is not created
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    # setup
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    dset.assign_groups({
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7]
    })
    add_log2fc(dset, ['Par', 'Dap2'])
    add_2group_pvalue(dset, ['Par', 'Dap2'], 'students')

    img_dir = os.path.dirname(__file__)
    fig_path = os.path.join(img_dir, 'volcano_Par-Dap2_students_raw.png')
    volcano_2group(dset, ['Par', 'Dap2'], 'students', img_dir)

    if not os.path.isfile(fig_path):
        m = 'volcano_2group_real1: image file {} not found'
        raise RuntimeError(m.format(fig_path))
    else:
        # delete the image file
        os.remove(fig_path)

    return True


def volcano_2group_badstats_real1():
    """
volcano_2group_badstats_real1
    description:
        Tests the function that makes a volcano plot by giving it bad information, expects ValueErrors

        Test fails if there are any errors other than the expected ValueErrors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    # setup
    dset = Dataset(os.path.join(os.path.dirname(__file__), 'real_data_1.csv'), esi_mode='neg')
    dset.assign_groups({
        'Par': [0, 1, 2, 3],
        'Dap2': [4, 5, 6, 7]
    })
    img_dir = os.path.dirname(__file__)
    fig_path = os.path.join(img_dir, 'volcano_Par-Dap2_students_raw.png')
    # first try using 3 group names
    try:
        volcano_2group(dset, ['Par', 'Dap2', 'Dal2'], 'students', img_dir)
    except ValueError:
        pass
    # now use the correct two groups without p-value calculated
    try:
        volcano_2group(dset, ['Par', 'Dap2'], 'students', img_dir)
    except ValueError:
        pass
    add_2group_pvalue(dset, ['Par', 'Dap2'], 'students')
    # now use the correct two groups and wrong stats test for p-value
    try:
        volcano_2group(dset, ['Par', 'Dap2'], 'welchs', img_dir)
    except ValueError:
        pass
    # now use the correct two groups with correct p-value and no log2fc
    try:
        volcano_2group(dset, ['Par', 'Dap2'], 'students', img_dir)
    except ValueError:
        pass

    return True


# references to al of the test functions to be run, and order to run them in
all_tests = [
    barplot_feature_bygroup_mock1,
    batch_barplot_feature_bygroup_real1,
    scatter_pca3_projections_bygroup_mock1,
    scatter_plsda_projections_bygroup_mock1,
    splot_plsda_pcorr_bygroup_mock1,
    scatter_plsra_projections_bygroup_real1,
    heatmap_lipid_class_log2fc_real1,
    volcano_2group_real1,
    volcano_2group_badstats_real1
]
if __name__ == '__main__':
    run_tests(all_tests)

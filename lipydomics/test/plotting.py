"""
    lipydomics/test/plotting.py
    Dylan H. Ross
    2019/09/25

    description:
        Define tests associated with the lipydomics/plotting.py module
"""


import os
import numpy as np

from lipydomics.data import Dataset
from lipydomics.stats import add_pca3, add_plsda, add_2group_corr
from lipydomics.plotting import (
    barplot_feature_bygroup, scatter_pca3_projections_bygroup, scatter_plsda_projections_bygroup, 
    splot_plsda_pcorr_bygroup
)


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
    img_path = os.path.join(img_dir, 'PLS-DA_projections_{}_normed.png'.format('-'.join(group_names)))

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
    img_path = os.path.join(img_dir, 'S-plot_{}_normed.png'.format('-'.join(group_names)))

    splot_plsda_pcorr_bygroup(dset, group_names, img_dir, normed=True)
    if not os.path.isfile(img_path):
        m = 'scatter_plsda_projections_bygroup_mock1: image file {} not found'
        raise RuntimeError(m.format(img_path))
    else:
        # delete the image file 
        os.remove(img_path)
        
    return True


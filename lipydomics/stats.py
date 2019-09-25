"""
    lipydomics/stats.py
    Dylan H. Ross
    2019/02/03

    description:
        A set of functions for performing statistical analyses on the lipidomics data. Generally, these functions 
        should produce one or more columns to associate with the data, as well as a label describing the analysis
        performed (what groupings are used, normalized or raw data, etc.). This data is added to the Dataset instance
        into a dictionary called stats, where the key is the label string and the value is a numpy.ndarray containing
        the desired statistic. Some functions will also return extra information, for example: when doing PCA analysis,
        the loadings of each feature are added as data in Dataset.stats and projections for each sample within user-
        specified groups can optionally be returned.

        if labels are not explicitly provided, they can generated in a systematic fashion:
            statistic_groups_raw-or-normed
        examples:
            ANOVA_ctrl-AY-VitE-AYVitE_normed
            PCA_ALL_normed

        WARNING: adding statistics with the same label as ones already added to the Dataset object will overwrite the 
        existing data
"""


import numpy as np
from scipy.stats import f_oneway
from sklearn.decomposition import PCA


def add_anova_p(dataset, group_names, normed=False):
    """
add_anova_p
    description:
        adds a column containing ANOVA p-values computed for user-specified groups
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- groups to use to compute the ANOVA p-value
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
"""
    group_data = np.array([_ for _ in dataset.get_data_bygroup(group_names, normed=normed)])
    anova_p = np.array([f_oneway(*feature)[1] for feature in zip(*group_data)])
    if normed:
        normed = "normed"
    else: 
        normed = "raw"
    label = "ANOVA_" + '-'.join(group_names) + "_" + normed
    # add the statistic into the Dataset
    dataset.stats[label] = anova_p


def add_pca3(dataset, normed=False, random_state=69):
    """
add_pca3
    description:
        Computes a 3-component PCA from all features in the Dataset (using either raw or normalized intensities) and
        adds the relevant information (feature loadings, projections) to the Dataset. 

        The fitted PCA object is added as an instance variable (Dataset.pca3_). The feature loadings (3, n_features) 
        and PCA projections (n_samples, 3) are added to Dataset.stats with the labels 'PCA3_loadings_{raw/normed}' and 
        'PCA3_projections_{raw/normed}', respectively.
        
        If the Dataset.pca3_ instance variable is already present, then it will be overridden by the new one generated 
        by this function.
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
        [random_state (int)] -- pRNG seed for deterministic results [optional, default=69]
"""
    # initialize the PCA object, add it to the Dataset
    dataset.pca3_ = PCA(n_components=3, random_state=random_state)
    # fit the PCA to the data
    if normed:
        nrm = 'normed'
        # need to transpose the intensities array to shape: (n_samples, n_features) for PCA
        use_data = dataset.normed_intensities.T
    else:
        nrm = 'raw'
        # need to transpose the intensities array to shape: (n_samples, n_features) for PCA
        use_data = dataset.intensities.T
    dataset.pca3_.fit(use_data)
    # add the statistic into the Dataset
    dataset.stats['PCA3_loadings_{}'.format(nrm)] = dataset.pca3_.components_
    dataset.stats['PCA3_projections_{}'.format(nrm)] = dataset.pca3_.transform(use_data)


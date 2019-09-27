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
from scipy.stats import f_oneway, pearsonr
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression


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
        
        If the Dataset.pca3_ instance variable or either of the Dataset.stats entries are already present, then they 
        will be overridden.
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

    # add the statistics into the Dataset
    dataset.stats['PCA3_loadings_{}'.format(nrm)] = dataset.pca3_.components_
    dataset.stats['PCA3_projections_{}'.format(nrm)] = dataset.pca3_.transform(use_data)


def add_plsda(dataset, group_names, normed=False):
    """
add_plsda
    description:
        Performs PLS-DA using two specified groups (e.g. groups A and B).

        Uses the PLSRegression class from sklearn with the target variable being simply 1 for group A or -1 for group B
        which effectively converts the task from regression to calssification. The fitted PLSRegression object is stored 
        in an instance variable (Dataset.plsr_). The feature loadings (n_features, 2) and projections (n_samples, 2) are 
        added to Dataset.stats with the labels 'PLS-DA_A-B_loadings_{raw/normed}' and 
        'PLS-DA_A-B_projections_{raw/normed}', respectively.
    
        If the Dataset.plsr_ instance variable or either of the Dataset.stats entries are already present, then they 
        will be overridden.
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- groups to use to compute the PLS-DA, only 2 groups allowed
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
"""
    if len(group_names) != 2:
        m = 'add_plsda: 2 group names must be specified for PLS-DA, {} group names specified'
        raise ValueError(m.format(len(group_names)))

    # get the group data, reshape and concatenate -> X 
    X = np.concatenate([_.T for _ in dataset.get_data_bygroup(group_names, normed=normed)])
    # target variable is just 1 for group A and -1 for group B
    n_A = len(dataset.group_indices[group_names[0]])
    n_B = len(dataset.group_indices[group_names[1]])
    y = np.array([1 for _ in range(n_A)] + [-1 for _ in range(n_B)])

    # initialize the PLSRegression object, add to the Dataset, fit the group data
    dataset.plsr_ = PLSRegression()
    dataset.plsr_.fit(X, y)

    if normed:
        nrm = 'normed'
    else:
        nrm = 'raw'

    # add the statistics into the Dataset
    dataset.stats['PLS-DA_{}_loadings_{}'.format('-'.join(group_names), nrm)] = dataset.plsr_.x_loadings_
    dataset.stats['PLS-DA_{}_projections_{}'.format('-'.join(group_names), nrm)] = dataset.plsr_.transform(X)


def add_2group_corr(dataset, group_names, normed=False):
    """
add_2group_corr
    description:
        Computes Pearson correlation coefficient between two specified groups (e.g. groups A and B) for all features.

        The Pearson correlation coefficients (n_features,) are added to Dataset.stats with the label 
        '2-group-corr_A-B_{raw/normed}'. These correlation coefficients can be combined with the loadings from a PLS-DA
        computed on the same groups in order to generate the familiar S-plot.
    
        * To use these data with the loadings from PLS-DA, the same group names must be specified in the same order *

        If the Dataset.stats entry is already present, then it will be overridden.
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- groups to use to compute the PLS-DA, only 2 groups allowed
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
"""
    if len(group_names) != 2:
        m = 'add_2group_corr: 2 group names must be specified for correlation, {} group names specified'
        raise ValueError(m.format(len(group_names)))

    # get the group data, reshape and concatenate -> X 
    Y = np.concatenate([_.T for _ in dataset.get_data_bygroup(group_names, normed=normed)]).T
    # target variable is just 1 for group A and -1 for group B
    n_A = len(dataset.group_indices[group_names[0]])
    n_B = len(dataset.group_indices[group_names[1]])
    x = np.array([1 for _ in range(n_A)] + [-1 for _ in range(n_B)])

    # compute correlation coefficients for each feature
    """
    corr = []
    for y in Y:
        c = 0.
        for _ in y:
            # if x is all 0s then dont bother computing the correlation coefficient, it just causes a warning
            if _ > 0:
                c = pearsonr(x, y)[0]
                break
        corr.append(c)

    # convert corr to a numpy array
    corr = np.array(corr)
    """
    corr = np.array([pearsonr(x, y)[0] for y in Y])

    if normed:
        nrm = 'normed'
    else:
        nrm = 'raw'

    # add the statistic into the Dataset
    dataset.stats['2-group-corr_{}_{}'.format('-'.join(group_names), nrm)] = corr
   

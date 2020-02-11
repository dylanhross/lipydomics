"""
    lipydomics/plotting.py
    Dylan H. Ross
    2019/09/23
    description:
        A set of functions for generating plots of various features, statistical analyses, etc.
"""


import os
import numpy as np
from csv import reader
from matplotlib import pyplot as plt
from matplotlib import rcParams


rcParams['font.size'] = 8
CS = ['#2FA2AB', '#9BD0B9', 'Purple', 'Blue', 'Green', 'Orange', 'Red', 'Yellow',
      '#E8ACF6', 'Grey', '#D6BF49', '#412F88', '#A2264B', '#3ACBE8', '#1CA3DE', '#0D85D8']
IMG_RES = 350  # image resolution


def barplot_feature_bygroup(dataset, group_names, img_dir, feature, normed=False, tolerance=(0.01, 0.1, 3.)):
    """
barplot_feature_bygroup
    description:
        generates a bar plot of the specified feature, comparing the mean intensities of the specified groups and saves
        the image to a specified directory. The filename of the image is:
            'bar_{mz}-{rt}-{ccs}_{group_name1}-{group_name2}-{etc.}_{raw or normed}.png'
        returns a boolean indicating whether a feature was found (and therefore a plot was generated)
        * m/z tolerance is in Da, rt tolerance is in min, CCS tolerance is in percent *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups to plot against
        img_dir (str) -- directory to save the image under
        feature (tuple(float)) -- m/z, rt, and ccs of feature
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
        [tolerance (tuple(float))] -- tolerance to use for m/z, rt, and ccs search, CCS tolerance is a percentage not
                                        an absolute tolerance [optional, default=(0.01, 0.1, 3.)]
    returns:
        (bool) -- at least one feature was found
"""
    # search for the corresponding feature
    mz_r, rt_r, ccs_r = feature
    mz_t, rt_t, ccs_t = tolerance
    found_feat = []
    for i in range(len(dataset.labels)):
        mz, rt, ccs = dataset.labels[i]
        if abs(mz - mz_r) <= mz_t and abs(rt - rt_r) <= rt_t and (100. * (abs(ccs - ccs_r) / ccs)) <= ccs_t:
            # matched a feature
            found_feat.append(i)

    if not found_feat:
        # we should not be printing messages at this level... let the interface handle that, this is why we have a
        # return status to check
        #print('did not find a feature matching mz: {:.4f} rt: {:.2f} ccs: {:.1f}'.format(mz_r, rt_r, ccs_r),
        #      'within tolerances:', tolerance)
        return False

    # go through each matched feature and generate a plot
    for i in found_feat:
        mz, rt, ccs = dataset.labels[i]
        # check for identifications first
        if dataset.feat_ids is not None:
            put_id, put_lvl = dataset.feat_ids[i], dataset.feat_id_levels[i]
        else:
            put_id, put_lvl = '', ''
        if type(put_id) == list:
            put_id = put_id[0] 

        # generate a filename
        if normed:
                nrm = 'normed'
        else: 
            nrm = 'raw'
        fig_name = 'bar_{:.4f}-{:.2f}-{:.1f}_'.format(mz, rt, ccs) + '-'.join(group_names) + '_{}.png'.format(nrm)
        fig_path = os.path.join(img_dir, fig_name)

        group_data = np.array([_[i] for _ in dataset.get_data_bygroup(group_names, normed=normed)])

        x = [_ for _ in range(len(group_data))]
        y = [np.mean(_) for _ in group_data]
        e = [np.std(_) for _ in group_data]
        c = [c_ for _, c_ in zip(x, CS)]

        fig = plt.figure(figsize=(1. + 0.5 * len(group_names), 2))
        ax = fig.add_subplot(111)

        ax.bar(x, y, yerr=e, color=c, width=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(group_names, rotation='vertical')
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
        ax.set_ylabel('intensity')
        
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

        ttl = 'mz: {:.4f} rt: {:.2f} ccs: {:.1f}'.format(mz, rt, ccs)
        ttl += '\n{} ({})'.format(put_id, put_lvl)
        ax.set_title(ttl, fontsize=6, y=1.075)

        plt.tight_layout()
        plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
        plt.close()
    # if we made it here, at least one feature was found
    return True


def batch_barplot_feature_bygroup(dataset, group_names, img_dir, in_csv, normed=False, tolerance=(0.01, 0.1, 3.)):
    """
barplot_feature_bygroup
    description:
        generates bar plots of features from an input .csv file, comparing the mean intensities of the specified groups
        and saves the images to a specified directory. The filename of the images are:
            'bar_{mz}-{rt}-{ccs}_{group_name1}-{group_name2}-{etc.}_{raw or normed}.png'
        * m/z tolerance is in Da, rt tolerance is in min, CCS tolerance is in percent *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups to plot against
        img_dir (str) -- directory to save the image under
        in_csv (str) -- filename of a .csv file with features to search for
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
        [tolerance (tuple(float))] -- tolerance to use for m/z, rt, and ccs search, CCS tolerance is a percentage not
                                        an absolute tolerance [optional, default=(0.01, 0.1, 3.)]
"""
    with open(in_csv, 'r') as inf:
        next(inf)  # expects  a header row
        rdr = reader(inf)
        for mz, rt, ccs in rdr:  # iterate through query values
            feat = (float(mz), float(rt), float(ccs))
            # make successive calls to barplot_feature_bygroup
            barplot_feature_bygroup(dataset, group_names, img_dir, feat, normed=normed, tolerance=tolerance)


def scatter_pca3_projections_bygroup(dataset, group_names, img_dir, normed=False):
    """
scatter_pca3_projections_bygroup
    description:
        generates a scatter plot of the PCA projections for a specified set of groups and saves the image to a 
        specified directory. The filename of the image is:
            'PCA3_{group_name1}-{group_name2}-{etc.}_projections_{raw or normed}.png'
        * The same group names (in the same order) as were used in the call to add_pca3(...) must be used. *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups used to calculate the PCA
        img_dir (str) -- directory to save the image under
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
"""
    # generate the path to save the figure under
    if normed:
        nrm = 'normed'
    else: 
        nrm = 'raw'
    fig_name = 'PCA3_{}_projections_{}.png'.format('-'.join(group_names), nrm)
    fig_path = os.path.join(img_dir, fig_name)

    # make the plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)

    ax.axvline(lw=0.5, c='k', zorder=0)
    ax.axhline(lw=0.5, c='k', zorder=0)

    d = dataset.stats['PCA3_{}_projections_{}'.format('-'.join(group_names), nrm)]
    si = []
    i = 0
    for gn in group_names:
        l = len(dataset.group_indices[gn])
        si.append(len(dataset.group_indices[gn]) + i)
        i += l

    for dg, c, gn in zip(np.split(d, si), CS, group_names):
        ax.scatter(*dg.T[:2], marker='.', s=24, c=c, label=gn)

    for d in ['top', 'right', 'bottom', 'left']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('PC1 ({:.1f} %)'.format(100. * dataset.pca3_.explained_variance_ratio_[0]), fontsize=8)
    ax.set_ylabel('PC2 ({:.1f} %)'.format(100. * dataset.pca3_.explained_variance_ratio_[1]), fontsize=8)
    ax.set_title('3 component PCA', fontsize=8, fontweight='bold')

    ax.ticklabel_format(style='sci', scilimits=(0, 0))

    ax.legend(fontsize=6, borderpad=0.2)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
    plt.close()


def scatter_plsda_projections_bygroup(dataset, group_names, img_dir, normed=False):
    """
scatter_plsda_projections_bygroup
    description:
        generates a scatter plot of the PLS-DA projections for a specified set of groups and saves the image to a 
        specified directory. The filename of the image is:
            'PLS-DA_projections_{group_A}-{group_B}_{raw or normed}.png'
        * The same group names (in the same order) as were used in the call to add_plsda(...) must be used. *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups to plot against
        img_dir (str) -- directory to save the image under
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
"""
    if len(group_names) != 2:
        m = 'scatter_plsda_projections_bygroup: 2 group names must be specified for PLS-DA, {} group names specified'
        raise ValueError(m.format(len(group_names)))

    # generate the path to save the figure under
    if normed:
        nrm = 'normed'
    else: 
        nrm = 'raw'
    fig_name = 'PLS-DA_projections_{}_{}.png'.format('-'.join(group_names), nrm)
    fig_path = os.path.join(img_dir, fig_name)

    # make the plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)

    ax.axvline(lw=0.5, c='k', zorder=0)
    ax.axhline(lw=0.5, c='k', zorder=0)

    #get the projections, split into groups A and B
    d = dataset.stats['PLS-DA_{}_projections_{}'.format('-'.join(group_names), nrm)]
    d_A = d[:len(dataset.group_indices[group_names[0]])]
    d_B = d[len(dataset.group_indices[group_names[0]]):]
    ax.scatter(*d_A.T, marker='.', s=24, c='r', label=group_names[0])
    ax.scatter(*d_B.T, marker='.', s=24, c='b', label=group_names[1])

    for d in ['top', 'right', 'bottom', 'left']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('scores[0]', fontsize=8)
    ax.set_ylabel('scores[1]', fontsize=8)
    ax.set_title('PLS-DA', fontsize=8, fontweight='bold')

    ax.ticklabel_format(style='sci', scilimits=(0, 0))

    ax.legend(fontsize=6, borderpad=0.2)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
    plt.close()


def splot_plsda_pcorr_bygroup(dataset, group_names, img_dir, normed=False):
    """
splot_plsda_pcorr_bygroup
    description:
        Generate an S-Plot using x_loadings from PLS-DA and Pearson correlation coefficients, and save the image to a 
        specified directory. The filename of the image is:
            'S-plot_{group_A}-{group_B}_{raw or normed}.png'
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups to plot against
        img_dir (str) -- directory to save the image under
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]   
"""
    if len(group_names) != 2:
        m = 'splot_plsda_pcorr_bygroup: 2 group names must be specified for S-plot, {} group names specified'
        raise ValueError(m.format(len(group_names)))

    # generate the path to save the figure under
    if normed:
        nrm = 'normed'
    else: 
        nrm = 'raw'
    fig_name = 'S-Plot_{}_{}.png'.format('-'.join(group_names), nrm)
    fig_path = os.path.join(img_dir, fig_name)

    # get the data
    x = dataset.stats['PLS-DA_{}_loadings_{}'.format('-'.join(group_names), nrm)].T[0]
    y = dataset.stats['2-group-corr_{}_{}'.format('-'.join(group_names), nrm)]

    # color the positive and negative values differently
    c = []
    for _y in y:
        if _y > 0:
            c.append('b')
        else:
            c.append('r')

    # make the plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)

    ax.axvline(lw=0.5, c='k', ls='--', zorder=0)
    ax.axhline(lw=0.5, c='k', ls='--', zorder=0)
    
    ax.scatter(x, y, c=c, s=1)

    ax.set_ylim([-1, 1])
    ax.set_xlabel('x loadings', fontsize=8)
    ax.set_ylabel('pcorr', fontsize=8)
    ax.set_title('S-plot\n(-1={}, 1={})'.format(*group_names), fontsize=8, fontweight='bold')
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))

    plt.tight_layout()
    plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
    plt.close()


def scatter_plsra_projections_bygroup(dataset, group_names, img_dir, normed=False):
    """
scatter_plsra_projections_bygroup
    description:
        generates a scatter plot of the PLS-RA projections for a specified set of groups and saves the image to a
        specified directory. The filename of the image is:
            'PLS-RA_{group_name1}-{group_name2}-{etc.}_projections_{raw or normed}.png'
        * The same group names (in the same order) as were used in the call to add_plsra(...) must be used. *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups used to calculate the PCA
        img_dir (str) -- directory to save the image under
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
"""
    # generate the path to save the figure under
    if normed:
        nrm = 'normed'
    else:
        nrm = 'raw'
    fig_name = 'PLS-RA_{}_projections_{}.png'.format('-'.join(group_names), nrm)
    fig_path = os.path.join(img_dir, fig_name)

    # make the plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)

    ax.axvline(lw=0.5, c='k', zorder=0)
    ax.axhline(lw=0.5, c='k', zorder=0)

    d = dataset.stats['PLS-RA_{}_projections_{}'.format('-'.join(group_names), nrm)]
    si = []
    i = 0
    for gn in group_names:
        l = len(dataset.group_indices[gn])
        si.append(len(dataset.group_indices[gn]) + i)
        i += l

    for dg, c, gn in zip(np.split(d, si), CS, group_names):
        ax.scatter(*dg.T[:2], marker='.', s=24, c=c, label=gn)

    for d in ['top', 'right', 'bottom', 'left']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('scores[0]', fontsize=8)
    ax.set_ylabel('scores[1]', fontsize=8)
    ax.set_title('PLS-RA', fontsize=8, fontweight='bold')

    ax.ticklabel_format(style='sci', scilimits=(0, 0))

    ax.legend(fontsize=6, borderpad=0.2)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
    plt.close()

"""
    lipydomics/plotting.py
    Dylan H. Ross
    2019/09/23
    description:
        A set of functions for generating plots of various features, statistical analyses, etc.
"""


import os
import numpy as np
from scipy.sparse import coo_matrix
from csv import reader
from matplotlib import pyplot as plt, colors as mcolors
from matplotlib import rcParams

from lipydomics.util import fetch_lipid_class_log2fc


rcParams['font.size'] = 8
CS = ['#2FA2AB', '#9BD0B9', 'Purple', 'Blue', 'Green', 'Orange', 'Red', 'Yellow',
      '#E8ACF6', 'Grey', '#D6BF49', '#412F88', '#A2264B', '#3ACBE8', '#1CA3DE', '#0D85D8']
IMG_RES = 350  # image resolution


def barplot_feature_bygroup(dataset, group_names, img_dir, feature, normed=False, tolerance=(0.01, 0.1, 1.)):
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
                                        an absolute tolerance [optional, default=(0.01, 0.1, 1.)]
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
    fig = plt.figure(figsize=(2.5, 2.5))
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
    fig_name = 'PLS-DA_{}_projections_{}.png'.format('-'.join(group_names), nrm)
    fig_path = os.path.join(img_dir, fig_name)

    # make the plot
    fig = plt.figure(figsize=(2.5, 2.5))
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
    fig = plt.figure(figsize=(2.5, 2.5))
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


def heatmap_lipid_class_log2fc(lipid_class, dataset, group_names, img_dir, normed=False):
    """
heatmap_lipid_class_log2fc
    description:
        generates a heatmap of log2(fold-chage) for all identified lipids in a lipid class and saves the image to a
        specified directory. The filename of the image is:
            '{lipid_class}_{group_name1}-{group_name2}_{etc.}_log2fc_{raw or normed}.png'
        * The same group names (in the same order) as were used in the call to add_log2fc(...) must be used. *
    parameters:
        lipid_class (str) -- lipid class to analyze
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- groups used to compute fold-change, only 2 groups allowed
        img_dir (str) -- directory to save the image under
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
    returns:
        (bool) -- found data for the specified lipid class
"""
    # check that identifications have been made first
    if dataset.feat_ids is None:
        m = 'heatmap_lipid_class_log2fc: no lipid identifications have been made'
        raise ValueError(m)

    # check that log2fa has been computed
    log2fa_label = 'LOG2FC_{}_{}'.format('-'.join(group_names), 'normed' if normed else 'raw')
    if log2fa_label not in dataset.stats:
        m = 'heatmap_lipid_class_log2fc: required statistic "{}" not in dataset.stats'
        raise ValueError(m.format(log2fa_label))

    # assemble the data
    nc, nu, l2 = fetch_lipid_class_log2fc(lipid_class, dataset, group_names, normed=normed)
    if nc is None:
        # no matching data was found
        return False

    # get the min/max n_carbon, n_unsat, and log2fc
    min_nc, max_nc = np.min(nc), np.max(nc)
    min_nu, max_nu = np.min(nu), np.max(nu)
    max_abs_l2 = np.max(np.abs(l2))
    # create an appropriately sized sparse matrix to hold the data
    spmat = coo_matrix((l2, (nc, nu)))

    # generate the figure
    fig = plt.figure(figsize=(3.33, 3.33))
    ax = fig.add_subplot(111)
    # plot the fold-changes
    im = ax.pcolor(spmat.toarray(), cmap='bwr', vmin=-max_abs_l2, vmax=max_abs_l2, edgecolors='k', linewidth=0.5)
    fig.colorbar(im, label='log2({} / {})'.format(group_names[1], group_names[0]))
    # mask out the missing values
    mask_cm = mcolors.ListedColormap([np.array([0., 0., 0., 0.4]), np.array([0., 0., 0., 0.])], N=2)
    bnorm = mcolors.BoundaryNorm([0., 0.0001, 100000.], mask_cm.N, clip=True)
    ax.pcolor(np.abs(spmat.toarray()), cmap=mask_cm, norm=bnorm)
    # adjustments
    ax.set_ylim([min_nc, max_nc + 1])
    ax.set_xticks([_ for _ in range(0, max_nu + 1, 2)])
    ax.set_title(lipid_class)
    ax.set_xlabel('FA unsaturations')
    ax.set_ylabel('FA carbons')
    plt.tight_layout()
    # save the figure
    # generate the path to save the figure under
    fig_name = '{}_{}_log2fc_{}.png'.format(lipid_class, '-'.join(group_names), 'normed' if normed else 'raw')
    fig_path = os.path.join(img_dir, fig_name)
    plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
    plt.close()

    # everything worked
    return True


def volcano_2group(dataset, group_names, stats_test, img_dir, normed=False, p_upper_bound=0.01):
    """
heatmap_lipid_class_log2fc
    description:
        generates a volcano plot (Log2(FC) vs. -Log10(p-value)) for a two-group comparison and saves the image to a
        specified directory. The filename of the image is:
            'volcano_{group_name1}-{group_name2}_{stats_test}_{raw or normed}.png'
        * The same group names (in the same order) as were used in the call to add_log2fc(...) and the call to
          add_2group_pvalue(...) must be used. *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- groups to use to compute fold-change, only 2 groups allowed
        stats_test (str) -- specify the statistical test that was used
        img_dir (str) -- directory to save the image under
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
        [p_lower_bound (float)] -- upper limit for p-value, lower p-values get colors [optional, default=0.01]
    returns:
        (bool) -- found data for the specified lipid class
"""
    if len(group_names) != 2:
        m = 'volcano_2group: 2 group names must be specified for volcano plot, {} group names specified'
        raise ValueError(m.format(len(group_names)))

    # grab the pvalues, take the negative log10
    nrm = 'normed' if normed else 'raw'
    stat_abbrev = {'students': 'studentsP', 'welchs': 'welchsP', 'mann-whitney': 'mannwhitP'}[stats_test]
    pv_label = "{}_".format(stat_abbrev) + '-'.join(group_names) + "_" + nrm
    if pv_label not in dataset.stats:
        m = 'volcano_2group: p-value data ("{}") not present in Dataset.stats'
        raise ValueError(m.format(pv_label))
    pv = dataset.stats[pv_label]
    l10pv = np.log10(pv) * -1.

    # grab the log2(fc)
    log2fa_label = 'LOG2FC_{}_{}'.format('-'.join(group_names), nrm)
    if log2fa_label not in dataset.stats:
        m = 'volcano_2group: Log2(FC) data ("{}") not present in Dataset.stats'
        raise ValueError(m.format(log2fa_label))
    l2fc = dataset.stats[log2fa_label]

    # generate the path to save the figure under
    fig_name = 'volcano_{}_{}_{}.png'.format('-'.join(group_names), stats_test, nrm)
    fig_path = os.path.join(img_dir, fig_name)

    # make the plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)

    # color each of the points according to their values
    c = []
    p_ub = -1. * np.log10(p_upper_bound)
    for l2fc_, l10pv_ in zip(l2fc, l10pv):
        # any p-value greater than 0.05 is grey
        if l10pv_ < p_ub:
            c.append('grey')
        else:
            if l2fc_ > 0:
                # increased features are red
                c.append('r')
            else:
                # decreased features are blue
                c.append('b')

    ax.axvline(lw=0.5, ls='--', c='k', zorder=0)
    ax.axhline(p_ub, lw=0.5, ls='--', c='grey', zorder=0)

    ax.scatter(l2fc, l10pv, s=4, c=c, edgecolors='none')

    ax.set_xlabel(r'$Log_2(FC)$', fontsize=8)
    ax.set_ylabel(r'$-Log_{10}(Pvalue)$', fontsize=8)
    ax.set_title('{} vs. {}'.format(*group_names), fontsize=8, fontweight='bold')
    plt.tight_layout()
    plt.savefig(fig_path, dpi=IMG_RES, bbox_inches='tight')
    plt.close()



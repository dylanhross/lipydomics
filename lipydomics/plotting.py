"""
    lipydomics/plotting.py
    Dylan H. Ross
    2019/09/23

    description:
        A set of functions for generating plots of various features, statistical analyses, etc.

"""


import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams


rcParams['font.size'] = 8


def barplot_feature_bygroup(dataset, group_names, path, feature, normed=False, tolerance=(0.01, 0.1, 1.)):
    """
barplot_feature_bygroup
    description:
        generates a bar plot of the specified feature, comparing the mean intensities of the specified groups and saves
        the image to a specified directory. The filename of the image is:
            'bar_{mz}-{rt}-{ccs}_{group_name1}-{group_name2}-{etc.}_{raw or normed}.png'
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- pick groups to plot against
        path (str) -- path to save the image under
        feature (tuple(float)) -- m/z, rt, and ccs of feature
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
        [tolerance (tuple(float))] -- tolerance to use for m/z, rt, and ccs search [optional, default=(0.01, 0.1, 1.)] 
"""
    # search for the corresponding feature
    mz_r, rt_r, ccs_r = feature
    mz_t, rt_t, ccs_t = tolerance
    found_feat = []
    for i in range(len(dataset.labels)):
        mz, rt, ccs = dataset.labels[i]
        if abs(mz - mz_r) <= mz_t and abs(rt - rt_r) <= rt_t and abs(ccs - ccs_r) <= ccs_t:
            # matched a feature
            found_feat.append(i)

    if not found_feat:
        print('did not find a feature matching mz: {:.4f} rt: {:.2f} ccs: {:.1f}'.format(mz_r, rt_r, ccs_r), 
              'within tolerances:', tolerance)
        return None

    # go through each matched feature and generate a plot
    for i in found_feat:
        mz, rt, ccs = dataset.labels[i]
        # generate a filename
        if normed:
                nrm = 'normed'
        else: 
            nrm = 'raw'
        fig_name = 'bar_{:.4f}-{:.2f}-{:.1f}_'.format(mz, rt, ccs) + '-'.join(group_names) + '_{}.png'.format(nrm)
        fig_path = os.path.join(path, fig_name)

        group_data = np.array([_[i] for _ in dataset.get_data_bygroup(group_names, normed=normed)])

        x = [_ for _ in range(len(group_data))]
        y = [np.mean(_) for _ in group_data]
        e = [np.std(_) for _ in group_data]
        c = [c_ for _, c_ in zip(x, ['r', 'b', '#ffa600', 'purple', 'green', 'm'])]

        fig = plt.figure(figsize=(1. + 0.5 * len(group_names), 2))
        ax = fig.add_subplot(111)

        ax.bar(x, y, yerr=e, color=c, width=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(group_names, rotation='vertical')
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
        ax.set_ylabel('intensity')
        ax.set_title('mz: {:.4f} rt: {:.2f} ccs: {:.1f}'.format(mz, rt, ccs), fontsize=8)

        plt.tight_layout()
        plt.savefig(fig_path, dpi=200, bbox_inches='tight')



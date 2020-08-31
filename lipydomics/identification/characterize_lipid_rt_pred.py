"""
    lipydomics/identification/characterize_lipid_rt_pred.py
    Dylan H. Ross
    2019/12/03

    description:
        Characterizes performance of the predictive model for generating predicted RT values by generating plots
        of predicted vs. measured RT
"""


from sqlite3 import connect
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from numpy import median

from ..util import print_and_log
from .encoder_params import rt_lipid_classes, rt_fa_mods


rcParams['font.size'] = 6


def single_class_plot(cursor, lipid_class, fa_mod=None):
    """
single_class_plot
    description:
        generates a plot comparing predicted RT values against experimentally measured ones for a given lipid class
        (accounting for fa_mod if present). Saves the plot as a figure in the rt_pred_perf/ directory
    parameters:
        cursor (sqlite3.cursor) -- cursor for querying lipids.db
        lipid_class (str) -- lipid class
        [fa_mod (None or str)] -- fatty acid modifier [optional, default=None]
"""
    # predicted and measured retention time
    rt_t, rt_m = [], []

    # fetch the predicted data
    if fa_mod in rt_fa_mods:
        qry = 'SELECT rt FROM predicted_mz JOIN predicted_rt ON predicted_mz.t_id=predicted_rt.t_id '
        qry += 'WHERE lipid_class="{}" AND fa_mod=="{}"'
        qry = qry.format(lipid_class, fa_mod)
    else:
        qry = 'SELECT rt FROM predicted_mz JOIN predicted_rt ON predicted_mz.t_id=predicted_rt.t_id '
        qry += 'WHERE lipid_class="{}" AND fa_mod IS NULL'
        qry = qry.format(lipid_class)
    for rt in cursor.execute(qry).fetchall():
        rt_t.append(float(rt[0]))

    # fetch the measured data
    if fa_mod in rt_fa_mods:
        qry = 'SELECT rt FROM measured WHERE lipid_class="{}" AND fa_mod=="{}" AND rt IS NOT NULL'
        qry = qry.format(lipid_class, fa_mod)
    else:
        qry = 'SELECT rt FROM measured WHERE lipid_class="{}" AND fa_mod IS NULL AND rt IS NOT NULL'
        qry = qry.format(lipid_class)
    for rt in cursor.execute(qry).fetchall():
        rt_m.append(float(rt[0]))

    this_dir = os.path.dirname(__file__)
    fig_fname = 'rt_pred_perf/{}_{}{}.png'.format(len(rt_m), lipid_class, fa_mod if fa_mod else '')
    fig_path = os.path.join(this_dir, fig_fname)

    if len(rt_m) > 0 and len(rt_t) > 0:

        # determine the median of all of the RT values and set the y-bounds of the plots accordingly
        rt_all = rt_m + rt_t
        rt_min, rt_max = min(rt_all), max(rt_all)
        rt_middle = round((rt_max - rt_min) / 2. + rt_min, 1)

        fig = plt.figure(figsize=(1.2, 1.8))
        ax = fig.add_subplot(111)

        bp1 = ax.boxplot(rt_t, positions=[1], widths=0.8)
        bp2 = ax.boxplot(rt_m, positions=[2], widths=0.8, flierprops={'markersize': 3})

        plt.setp(bp1['whiskers'], color='r')
        plt.setp(bp2['whiskers'], color='b')
        plt.setp(bp1['boxes'], color='r')
        plt.setp(bp2['boxes'], color='b')
        plt.setp(bp1['caps'], color='r')
        plt.setp(bp2['caps'], color='b')
        plt.setp(bp1['medians'], color='r')
        plt.setp(bp2['medians'], color='b')

        ax.set_xlim([0, 3])
        ax.set_xticks([1, 2])
        ax.set_ylim([rt_middle - 0.7, rt_middle + 0.7])
        ax.set_xticklabels(['pred', 'meas (n={})'.format(len(rt_m))], rotation=45)
        ax.set_ylabel('retention time (min)')
        ax.set_title('{}{}'.format(lipid_class, fa_mod if fa_mod else ''))

        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)

        plt.savefig(fig_path, dpi=400, bbox_inches='tight')
        plt.close()


def main(tstamp):
    """ main build function """

    # connect to database
    db_path = os.path.join(os.path.dirname(__file__), 'lipids.db')
    con = connect(db_path)
    cur = con.cursor()

    build_log = os.path.join(os.path.dirname(__file__), 'builds/build_log_{}.txt'.format(tstamp))
    with open(build_log, 'a') as bl:

        print_and_log('characterizing RT prediction performance ...', bl, end=' ')

        # automatically generate plots for all combinations
        qry = 'SELECT lipid_class, fa_mod FROM measured '
        qry += 'WHERE rt IS NOT NULL GROUP BY lipid_class, fa_mod'
        for lipid_class, fa_mod in cur.execute(qry).fetchall():
            # only use the classes and fa_mods that are explicitly encoded
            lc_ok = lipid_class in rt_lipid_classes
            fam_ok = fa_mod is None or fa_mod in rt_fa_mods
            if lc_ok and fam_ok:
                single_class_plot(cur, lipid_class, fa_mod=fa_mod)

        print_and_log('ok\n', bl)

    # close database connection
    con.close()


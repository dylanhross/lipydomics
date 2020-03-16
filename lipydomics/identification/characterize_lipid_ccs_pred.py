"""
    lipydomics/identification/characterize_lipid_ccs_pred.py
    Dylan H. Ross
    2019/10/09

    description:
        Characterizes performance of the predictive model for generating theoretical CCS values by generating plots
        of predicted vs. measured CCS organized by lipid class (along with FA modifier) and MS adduct
"""


import os
from sqlite3 import connect
from matplotlib import pyplot as plt
from matplotlib import rcParams

from .build_params import ccs_pred_ref_dsets
from .encoder_params import ccs_lipid_classes, ccs_fa_mods, ccs_ms_adducts


rcParams['font.size'] = 6


def single_class_plot(cursor, lipid_class, adduct, fa_mod=None):
    """
single_class_plot
    description:
        generates a plot comparing theoretical CCS values against experimentally measured ones for a given lipid class
        (accounting for fa_mod if present) and MS adduct. Saves the plot as a figure in the ccs_pred_perf/ directory
    parameters:
        cursor (sqlite3.cursor) -- cursor for querying lipids.db
        lipid_class (str) -- lipid class
        adduct (str) -- MS adduct
        [fa_mod (None or str)] -- fatty acid modifier [optional, default=None]
"""
    # theoretical and measured m/z and ccs
    mz_t, ccs_t = [], []
    mz_m, ccs_m = [], []

    # fetch the measured data
    if fa_mod in ['o', 'p']:
        qry = 'SELECT mz, ccs, src_tag FROM measured WHERE lipid_class="{}" AND fa_mod=="{}" AND adduct="{}"'
        qry = qry.format(lipid_class, fa_mod, adduct)
    else:
        qry = 'SELECT mz, ccs, src_tag FROM measured WHERE lipid_class="{}" AND fa_mod IS NULL AND adduct="{}"'
        qry = qry.format(lipid_class, adduct)
    for mz, ccs, st in cursor.execute(qry).fetchall():
        # only include data from the designated sources
        src_ok = st in ccs_pred_ref_dsets
        if src_ok:
            mz_m.append(float(mz))
            ccs_m.append(float(ccs))

    # if no measured data was found, return None to skip plotting this class
    if len(mz_m) < 1:
        return None

    # set bounds on the theoretical data to fetch and display
    mz_min, mz_max = int(min(mz_m)), int(max(mz_m))
    ccs_min, ccs_max = int(min(ccs_m)), int(max(ccs_m))

    # fetch the theoretical data
    if fa_mod in ['o', 'p']:
        qry = 'SELECT mz, ccs FROM theoretical_mz JOIN theoretical_ccs ON '
        qry += 'theoretical_mz.t_id=theoretical_ccs.t_id '
        qry += 'WHERE lipid_class="{}" AND fa_mod=="{}" AND adduct="{}" AND (mz BETWEEN {} AND {}) AND '
        qry += '(ccs BETWEEN {} AND {})'
        qry = qry.format(lipid_class, fa_mod, adduct, mz_min, mz_max, ccs_min, ccs_max)
    else:
        qry = 'SELECT mz, ccs FROM theoretical_mz JOIN theoretical_ccs ON '
        qry += 'theoretical_mz.t_id=theoretical_ccs.t_id '
        qry += 'WHERE lipid_class="{}" AND fa_mod IS NULL AND adduct="{}" AND (mz BETWEEN {} AND {}) AND '
        qry += '(ccs BETWEEN {} AND {})'
        qry = qry.format(lipid_class, adduct, mz_min, mz_max, ccs_min, ccs_max)
    for mz, ccs in cursor.execute(qry).fetchall():
        mz_t.append(float(mz))
        ccs_t.append(float(ccs))

    this_dir = os.path.dirname(__file__)
    fig_fname = 'ccs_pred_perf/{}_{}{}_{}.png'.format(len(ccs_m), lipid_class, fa_mod if fa_mod else '', adduct)
    fig_path = os.path.join(this_dir, fig_fname)

    fig = plt.figure(figsize=(3.33, 2))
    ax = fig.add_subplot(111)

    ax.scatter(mz_t, ccs_t, marker='.', s=8, c='#ffa600', label='theoretical')
    ax.scatter(mz_m, ccs_m, marker='.', s=1, c='purple', label='measured\n(n={} in training data)'.format(len(mz_m)))

    ax.legend()
    ax.set_xlabel('m/z')
    ax.set_ylabel(r'CCS ($\AA^2$)')
    ax.set_title('{}{} {}'.format(lipid_class, fa_mod if fa_mod else '', adduct), fontweight='bold')

    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)

    #ax.set_xlim([mz_min - 10, mz_max + 10])
    #ax.set_ylim([ccs_min - 2, ccs_max + 2])

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

        print('characterizing CCS prediction performance ...', end=' ')
        print('characterizing CCS prediction performance ...', end=' ', file=bl)

        # automatically generate plots for all data included in the model training
        qry = 'SELECT lipid_class, fa_mod, adduct, COUNT(*) as c FROM measured '
        qry += 'GROUP BY lipid_class, fa_mod, adduct HAVING c > 9'
        for lipid_class, fa_mod, adduct, c in cur.execute(qry).fetchall():
            # only use the classes, fa_mods and adducts that are explicitly encoded
            lc_ok = lipid_class in ccs_lipid_classes
            fam_ok = fa_mod is None or fa_mod in ccs_fa_mods
            add_ok = adduct in ccs_ms_adducts
            if lc_ok and fam_ok and add_ok:
                single_class_plot(cur, lipid_class, adduct, fa_mod=fa_mod)

        print('ok\n')
        print('ok\n', file=bl)

        # close database connection
        con.close()


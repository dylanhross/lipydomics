#!/usr/bin/python3
"""
    lipydomics/identification/characterize_lipid_ccs_pred.py
    Dylan H. Ross
    2019/10/09

    description:
        Characterizes performance of the predictive model for generating theoretical CCS values by generating plots
        of predicted vs. measured CCS organized by lipid class (along with FA modifier) and MS adduct
"""


from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 10


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
    fig_path = 'ccs_pred_perf/{}{}_{}.png'.format(lipid_class, fa_mod if fa_mod else '', adduct)

    # theoretical and measured m/z and ccs
    mz_t, ccs_t = [], []
    mz_m, ccs_m = [], []

    # fetch the theoretical data
    if fa_mod in ['o', 'p']:
        qry = 'SELECT mz, ccs FROM theoretical_mz JOIN theoretical_ccs ON theoretical_mz.t_id=theoretical_ccs.t_id '
        qry += 'WHERE lipid_class="{}" AND fa_mod=="{}" AND adduct="{}"'
        qry = qry.format(lipid_class, fa_mod, adduct)
    else:
        qry = 'SELECT mz, ccs FROM theoretical_mz JOIN theoretical_ccs ON theoretical_mz.t_id=theoretical_ccs.t_id '
        qry += 'WHERE lipid_class="{}" AND fa_mod IS NULL AND adduct="{}"'
        qry = qry.format(lipid_class, adduct)
    for mz, ccs in cursor.execute(qry).fetchall():
        mz_t.append(float(mz))
        ccs_t.append(float(ccs))

    # fetch the measured data
    if fa_mod in ['o', 'p']:
        qry = 'SELECT mz, ccs FROM measured WHERE lipid_class="{}" AND fa_mod=="{}" AND adduct="{}"'
        qry = qry.format(lipid_class, fa_mod, adduct)
    else:
        qry = 'SELECT mz, ccs FROM measured WHERE lipid_class="{}" AND fa_mod IS NULL AND adduct="{}"'
        qry = qry.format(lipid_class, adduct)
    for mz, ccs in cursor.execute(qry).fetchall():
        mz_m.append(float(mz))
        ccs_m.append(float(ccs))

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    ax.scatter(mz_t, ccs_t, marker='.', s=64, c='r', label='theo')
    ax.scatter(mz_m, ccs_m, marker='.', s=16, c='b', label='meas (n={})'.format(len(mz_m)))

    ax.legend()
    ax.set_xlabel('m/z')
    ax.set_ylabel('CCS (Å^2)')
    ax.set_title('{}{} {}'.format(lipid_class, fa_mod if fa_mod else '', adduct))

    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':

    from sqlite3 import connect
    import os

    # connect to the database
    con = connect('lipids.db')
    cur = con.cursor()

    print('characterizing CCS prediction performance ...', end=' ')

    # automatically generate plots for all combinations having at least 10 measured values
    qry = 'SELECT lipid_class, fa_mod, adduct, COUNT(*) as c FROM measured '
    qry += 'GROUP BY lipid_class, fa_mod, adduct HAVING c > 9'
    for lipid_class, fa_mod, adduct, c in cur.execute(qry).fetchall():
        single_class_plot(cur, lipid_class, adduct, fa_mod=fa_mod)

    print('ok\n')

    # close database connection
    con.close()


"""
    lipydomics/identification/id_levels.py
    Dylan H. Ross
    2020/03/31

    description:
        definitions for individual ID functions at various confidence levels
"""


from lipydomics.util import get_score


def id_feat_theo_mz(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm=None):
    """
id_feat_theo_mz
    description:
        identifies a feature on the basis of theoretical m/z
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (None)] -- NOT USED just included to maintain the same interface as the other id functions
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz FROM theoretical_mz WHERE mz BETWEEN ? AND ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz

    putative_ids, putative_scores = [], []
    for name, adduct, mz_x in cursor.execute(qry, (mz_min, mz_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs, mz_q=mz, mz_x=mz_x))

    if putative_ids:
        return putative_ids, 'theo_mz', putative_scores
    else:
        return '', '', []


def id_feat_meas_mz(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm=None):
    """
id_feat_meas_mz
    description:
        identifies a feature on the basis of measured m/z
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (None)] -- NOT USED just included to maintain the same interface as the other id functions
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz FROM measured WHERE mz BETWEEN ? AND ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz

    putative_ids, putative_scores = [], []
    for name, adduct, mz_x in cursor.execute(qry, (mz_min, mz_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs, mz_q=mz, mz_x=mz_x))

    if putative_ids:
        return putative_ids, 'theo_mz', putative_scores
    else:
        return '', '', []


def id_feat_theo_mz_rt(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
    """
id_feat_theo_mz_rt
    description:
        identifies a feature on the basis of theoretical m/z and retention time
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz, rt FROM theoretical_mz JOIN theoretical_rt ON ' \
            + 'theoretical_mz.t_id=theoretical_rt.t_id WHERE mz BETWEEN ? AND ? AND rt BETWEEN ? and ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    rt_min = rt - tol_rt
    rt_max = rt + tol_rt

    putative_ids, putative_scores = [], []
    for name, adduct, mz_x, rt_x in cursor.execute(qry, (mz_min, mz_max, rt_min, rt_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs, mz_q=mz, rt_q=rt, mz_x=mz_x, rt_x=rt_x))

    if putative_ids:
        return putative_ids, 'theo_mz_rt', putative_scores
    else:
        return '', '', []


def id_feat_theo_mz_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
    """
id_feat_theo_mz_ccs
    description:
        identifies a feature on the basis of theoretical m/z and CCS
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz, ccs FROM theoretical_mz JOIN theoretical_ccs ON ' \
            + 'theoretical_mz.t_id=theoretical_ccs.t_id WHERE mz BETWEEN ? AND ? AND ccs BETWEEN ? and ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    ccs_min = ccs - tol_ccs
    ccs_max = ccs + tol_ccs

    putative_ids, putative_scores = [], []
    for name, adduct, mz_x, ccs_x in cursor.execute(qry, (mz_min, mz_max, ccs_min, ccs_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs, mz_q=mz, ccs_q=ccs, mz_x=mz_x, ccs_x=ccs_x))

    if putative_ids:
        return putative_ids, 'theo_mz_ccs', putative_scores
    else:
        return '', '', []


def id_feat_theo_mz_rt_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
    """
id_feat_theo_mz_rt_ccs
    description:
        identifies a feature on the basis of theoretical m/z, retention time, and CCS
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz, rt, ccs FROM theoretical_mz JOIN theoretical_ccs ON ' \
            + 'theoretical_mz.t_id=theoretical_ccs.t_id JOIN theoretical_rt ON ' \
            + 'theoretical_mz.t_id=theoretical_rt.t_id WHERE mz BETWEEN ? AND ? AND ccs BETWEEN ? AND ? AND ' \
            + 'rt BETWEEN ? AND ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    ccs_min = ccs - tol_ccs
    ccs_max = ccs + tol_ccs
    rt_min = rt - tol_rt
    rt_max = rt + tol_rt

    putative_ids, putative_scores = [], []
    qdata = (mz_min, mz_max, ccs_min, ccs_max, rt_min, rt_max)
    for name, adduct, mz_x, rt_x, ccs_x in cursor.execute(qry, qdata).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs,
                                         mz_q=mz, rt_q=rt, ccs_q=ccs,
                                         mz_x=mz_x, rt_x=rt_x, ccs_x=ccs_x))

    if putative_ids:
        return putative_ids, 'theo_mz_rt_ccs', putative_scores
    else:
        return '', '', []


def id_feat_meas_mz_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
    """
id_feat_theo_meas_mz_ccs
    description:
        identifies a feature by matching a reference value on m/z, and CCS within tolerances
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz, ccs FROM measured WHERE mz BETWEEN ? AND ? AND ccs BETWEEN ? and ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    ccs_min = ccs - tol_ccs
    ccs_max = ccs + tol_ccs

    putative_ids, putative_scores = [], []
    qdata = (mz_min, mz_max, ccs_min, ccs_max)
    for name, adduct, mz_x, ccs_x in cursor.execute(qry, qdata).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs, mz_q=mz, ccs_q=ccs, mz_x=mz_x, ccs_x=ccs_x))

    if putative_ids:
        return putative_ids, 'meas_mz_ccs', putative_scores
    else:
        return '', '', []


def id_feat_meas_mz_rt(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
    """
id_feat_meas_mz_rt
    description:
        identifies a feature on the basis of measured m/z and retention time
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz, rt FROM measured WHERE mz BETWEEN ? AND ? AND rt BETWEEN ? and ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    rt_min = rt - tol_rt
    rt_max = rt + tol_rt

    putative_ids, putative_scores = [], []
    for name, adduct, mz_x, rt_x in cursor.execute(qry, (mz_min, mz_max, rt_min, rt_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs, mz_q=mz, rt_q=rt, mz_x=mz_x, rt_x=rt_x))

    if putative_ids:
        return putative_ids, 'theo_mz_rt', putative_scores
    else:
        return '', '', []


def id_feat_meas_mz_rt_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
    """
id_feat_theo_meas_mz_rt_ccs
    description:
        identifies a feature by matching a reference value on m/z, rt, and CCS within tolerances
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct, mz, rt, ccs FROM measured WHERE ' + \
            'mz BETWEEN ? AND ? AND rt BETWEEN ? and ? AND ccs BETWEEN ? and ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    rt_min = rt - tol_rt
    rt_max = rt + tol_rt
    ccs_min = ccs - tol_ccs
    ccs_max = ccs + tol_ccs

    putative_ids, putative_scores = [], []
    qdata = (mz_min, mz_max, rt_min, rt_max, ccs_min, ccs_max)
    for name, adduct, mz_x, rt_x, ccs_x in cursor.execute(qry, qdata).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))
        putative_scores.append(get_score(tol_mz, tol_rt, tol_ccs,
                                         mz_q=mz, rt_q=rt, ccs_q=ccs,
                                         mz_x=mz_x, rt_x=rt_x, ccs_x=ccs_x))

    if putative_ids:
        return putative_ids, 'meas_mz_rt_ccs', putative_scores
    else:
        return '', '', []


def id_feat_any(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2', use_rt=True):
    """
id_feat_any
    description:
        Goes through all defined levels of identification, starting with the highest level, then tries lower levels
        until an identification or identifications are made
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        mz (float) -- m/z to match
        rt (float) -- retention time to match
        ccs (float) -- CCS to match
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
        [use_rt (bool)] -- whether to use identification levels that involve retention time [optional, default=True]
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    if use_rt:
        id_funcs = [
            id_feat_meas_mz_rt_ccs,
            id_feat_theo_mz_rt_ccs,
            id_feat_meas_mz_rt,
            id_feat_theo_mz_rt,
            id_feat_meas_mz_ccs,
            id_feat_theo_mz_ccs,
            id_feat_meas_mz,
            id_feat_theo_mz
        ]
    else:
        id_funcs = [
            id_feat_meas_mz_ccs,
            id_feat_theo_mz_ccs,
            id_feat_meas_mz,
            id_feat_theo_mz
        ]

    for f in id_funcs:
        fid, lvl, scr = f(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm=norm)
        if fid:
            return fid, lvl, scr

    return '', '', []
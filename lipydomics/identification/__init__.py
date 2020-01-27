"""
    lipydomics/identification/__init__.py
    Dylan H. Ross
    2019/10/04

    description:
        Module for performing identification of individual lipid features
"""


from sqlite3 import connect
import os
from numpy import sqrt


def get_score(tol_mz, tol_rt, tol_ccs, mz_q=None, rt_q=None, ccs_q=None, mz_x=None, rt_x=None, ccs_x=None, norm='l2'):
    """
get_score
    description:
        computes a score reflecting the quality of an identification, using mz rt and ccs or any combination

        The score is determined by the residuals between the query values (q) and a potential match (x), normalized by
        their respective tolerances. If only a single pair of values (q and x) is provided, the score is simply the
        inverse of the normalized residual, otherwise, it is the inverse of the l1 or l2 norm of the normalized 
        residuals vector. The norm kwarg controls whether the l1 or l2 norm is used in computing the score 
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        tol_mz (float) -- tolerance for m/z 
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        [mz_q (None or float)] -- if specified, the query m/z [optional, default=None]
        [rt_q (None or float)] -- if specified, the query retention time [optional, default=None]
        [ccs_q (None or float)] -- if specified, the query CCS [optional, default=None]
        [mz_q (None or float)] -- if specified, the m/z of a potential match [optional, default=None]
        [rt_q (None or float)] -- if specified, the retention time of a potential match [optional, default=None]
        [ccs_q (None or float)] -- if specified, the CCS of a potential match [optional, default=None]
        [norm (str)] -- specify l1 or l2 norm [optional, default='l2']
    returns:
        (float) -- score (higher = more confidence in ID)
"""
    q = [mz_q, rt_q, ccs_q]
    x = [mz_x, rt_x, ccs_x]
    tol = [tol_mz, tol_rt, tol_ccs]
    # compute the normalized residuals
    rn = []
    for q_, x_, tol_ in zip(q, x, tol):
        if q_ and x_:
            rn.append((x_ - q_) / tol_)
    # if there was only a single value, just return the inverse of the absolute value
    if not rn:
        m = 'get_score: unable to compute residuals'
        raise RuntimeError(m)
    elif len(rn) == 1:
        # prevent zero-division just in case ...
        return 1. / max(abs(rn[0]), 0.000001)
    else:
        if norm == 'l1':
            return 1. / sum([abs(_) for _ in rn])
        elif norm == 'l2':
            return 1. / sqrt(sum([_**2. for _ in rn]))
        else:
            m = 'get_score: norm method "{}" not recognized'
            raise ValueError(m.format(norm))


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


def id_feat_any(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2'):
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
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    id_funcs = [
        id_feat_meas_mz_rt_ccs,
        id_feat_meas_mz_ccs,
        id_feat_theo_mz_rt_ccs,
        id_feat_theo_mz_ccs,
        id_feat_theo_mz_rt,
        id_feat_theo_mz
    ]

    for f in id_funcs:
        fid, lvl, scr = f(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode, norm='l2')
        if fid:
            return fid, lvl, scr

    return '', '', []


def add_feature_ids(dataset, tol, level='any', norm='l2'):
    """
add_feature_ids
    description:
        Goes through the list of features (mz, rt, ccs) in Dataset.labels and provides an identification for each.
        
        If a feature is unable to be identified using the specified criteria, then it is given a generic fixed-width 
        feature name based on mz, rt, and CCS:
            'UNK_{:09.4f}_{:05.2f}_{:06.2f}'.format(mz, rt, ccs)
            and the associated id_level is ''
        
        The method for making identifications is specified by the `level` param:
            (low)
            'theo_mz' -- simple matching based on theoretical m/z
            'theo_mz_ccs' -- match on theoretical m/z and CCS
            'theo_mz_rt_ccs' -- match on theoretical m/z, rt, and CCS
            'meas_mz_ccs' -- match on measured m/z and CCS
            'meas_mz_rt_ccs' -- match on m/z, rt, and CCS
            (high)
            'any' -- start at the highest level, then work downward

        The identifications are stored in `Dataset.feat_ids` and the associated identification levels are stored in 
        `Dataset.feat_id_levels`

        For each feature, the identification is a list of strings with all possible identifications (e.g. 
        ['PC(40:3)_[M+H]+', 'PE(42:2)_[M+H]+', 'PC(p36:1)_[M+Na]+']). The corresponding identification level is always 
        a single string (e.g. 'meas_mass_rt_ccs'). The identification score is a list of floats corresponding to scores
        for each of the possible identifications. The list is empty in the case of a failure to identify a feature.

        * whenever multiple identifications are returned, always sort in order of descending scores *

        The same tolerances are used for all levels of identification.

        Potential identifications are automatically restricted on the basis of electrospray ionization mode (taken from
        Dataset.esi_mode) if it is set to 'pos' or 'neg'

        > If a retention time calibration is available, use calibrated retention time for feature identification <

        * if identifications have already been made, subsequent calls to this function will override previous results *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        tol (tuple(float, float, float)) -- tolerance for m/z, rt, and CCS, respectively 
        [level (str)] -- specify the level of identification [optional, default='all']
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
"""
    if level not in ['theo_mz', 'theo_mz_ccs', 'theo_mz_rt_ccs', 'meas_mz_ccs', 'meas_mz_rt_ccs', 'any']:
        m = 'add_feature_ids: identification level "{}" not recognized'
        raise ValueError(m.format(level))

    # available identification functions
    id_funcs = {
        'any': id_feat_any,
        'meas_mz_rt_ccs': id_feat_meas_mz_rt_ccs,
        'meas_mz_ccs': id_feat_meas_mz_ccs,
        'theo_mz_rt_ccs': id_feat_theo_mz_rt_ccs,
        'theo_mz_ccs': id_feat_theo_mz_ccs,
        'theo_mz_rt': id_feat_theo_mz_rt,
        'theo_mz': id_feat_theo_mz
    } 

    # ESI mode from Dataset
    esi = dataset.esi_mode

    # initialize connection to lipids.db (stored within the lipydomics package)
    con = connect(os.path.join(os.path.dirname(__file__), 'lipids.db'))
    cur = con.cursor()

    feat_ids, feat_id_levels, feat_id_scores = [], [], []
    for mz, rt, ccs in dataset.labels:

        # use calibrated retention time if a retention time calibration has been set up
        rt = dataset.rt_calibration.get_calibrated_rt(rt) if dataset.rt_calibration is not None else rt

        # try to get identification(s)
        feat_id, feat_id_level, feat_id_score = id_funcs[level](cur, mz, rt, ccs, *tol, esi, norm=norm)

        if feat_id:
            if len (feat_id) > 1:
                # sort feat_id and feat_id_score in order of descending score
                feat_id, feat_id_score = [list(x) for x in zip(*sorted(zip(feat_id, feat_id_score), 
                                                                       key=lambda pair: pair[1], reverse=True))]
            feat_ids.append(feat_id)
            feat_id_levels.append(feat_id_level)
            feat_id_scores.append(feat_id_score)
        else:
            # fallback identification
            feat_id = 'UNK_{:09.4f}_{:05.2f}_{:06.2f}'.format(mz, rt, ccs)
            feat_ids.append(feat_id)
            feat_id_levels.append('')
            feat_id_scores.append([])

    # apply the identifications (and identification levels) to the Dataset
    dataset.feat_ids = feat_ids
    dataset.feat_id_levels = feat_id_levels
    dataset.feat_id_scores = feat_id_scores

    # close the database connection
    con.close()


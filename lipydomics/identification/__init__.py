"""
    lipydomics/identification/__init__.py
    Dylan H. Ross
    2019/10/04

    description:
        Module for performing identification of individual lipid features
"""


from sqlite3 import connect
import os


def id_feat_theo_mz(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode):
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
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct FROM theoretical_mz WHERE mz BETWEEN ? AND ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz

    putative_ids = []
    for name, adduct in cursor.execute(qry, (mz_min, mz_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))

    if putative_ids:
        return putative_ids, 'theo_mz'
    else:
        return '', ''


def id_feat_theo_mz_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode):
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
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct FROM theoretical_mz JOIN theoretical_ccs ON theoretical_mz.t_id=theoretical_ccs.t_id' \
            + ' WHERE mz BETWEEN ? AND ? AND ccs BETWEEN ? and ?'
    if esi_mode == 'pos':
        qry += ' AND adduct LIKE "%+"'
    elif esi_mode == 'neg':
        qry += ' AND adduct LIKE "%-"'

    mz_min = mz - tol_mz
    mz_max = mz + tol_mz
    ccs_min = ccs - tol_ccs
    ccs_max = ccs + tol_ccs

    putative_ids = []
    for name, adduct in cursor.execute(qry, (mz_min, mz_max, ccs_min, ccs_max)).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))

    if putative_ids:
        return putative_ids, 'theo_mz_ccs'
    else:
        return '', ''


def id_feat_meas_mz_rt_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode):
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
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    qry = 'SELECT name, adduct FROM measured WHERE mz BETWEEN ? AND ? AND rt BETWEEN ? and ? AND ccs BETWEEN ? and ?'
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

    putative_ids = []
    qdata = (mz_min, mz_max, rt_min, rt_max, ccs_min, ccs_max)
    for name, adduct in cursor.execute(qry, qdata).fetchall():
        putative_ids.append('{}_{}'.format(name, adduct))

    if putative_ids:
        return putative_ids, 'meas_mz_rt_ccs'
    else:
        return '', ''


def id_feat_any(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode):
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
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    return '', ''


def add_feature_ids(dataset, tol, level='any'):
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
            'meas_mz_rt_ccs' -- match on m/z, rt, and CCS for 
            (high)
            'any' -- start at the highest level, then work downward

        The identifications are stored in `Dataset.feat_ids` and the associated identification levels are stored in 
        `Dataset.feat_id_levels`

        For each feature, the identification is either a single string if there is only one identification (e.g. 
        'UNK_0123.4567_01.23_123.45') or a list of strings with all possible identifications (e.g. ['PC(40:3)_[M+H]+', 
        'PE(42:2)_[M+H]+', 'PC(p36:1)_[M+Na]+']). The corresponding identification level is always a single
        string (e.g. '' and 'meas_mass_rt_ccs', respectively, for the above examples) 

        The same tolerances are used for all levels of identification.

        Potential identifications are automatically restricted on the basis of electrospray ionization mode (taken from
        Dataset.esi_mode) if it is set to 'pos' or 'neg'

        * if identifications have already been made, subsequent calls to this function will override previous results *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        tol (tuple(float, float, float)) -- tolerance for m/z, rt, and CCS, respectively 
        [level (str)] -- specify the level of identification [optional, default='all']
"""
    if level not in ['theo_mz', 'theo_mz_ccs', 'meas_mz_rt_ccs', 'any']:
        m = 'add_feature_ids: identification level "{}" not recognized'
        raise ValueError(m.format(level))

    # available identification functions
    id_funcs = {
        'any': id_feat_any,
        'meas_mz_rt_ccs': id_feat_meas_mz_rt_ccs,
        'theo_mz_ccs': id_feat_theo_mz_ccs,
        'theo_mz': id_feat_theo_mz
    } 

    # ESI mode from Dataset
    esi = dataset.esi_mode

    # initialize connection to lipids.db (stored within the lipydomics package)
    con = connect(os.path.join(os.path.dirname(__file__), 'lipids.db'))
    cur = con.cursor()

    feat_ids, feat_id_levels = [], []
    for mz, rt, ccs in dataset.labels:

        # try to get identification(s)
        feat_id, feat_id_level = id_funcs[level](cur, mz, rt, ccs, *tol, esi)

        if feat_id:
            feat_ids.append(feat_id)
            feat_id_levels.append(feat_id_level)
        else:
            # fallback identification
            feat_id = 'UNK_{:09.4f}_{:05.2f}_{:06.2f}'.format(mz, rt, ccs)
            feat_ids.append(feat_id)
            feat_id_levels.append('')

    # apply the identifications (and identification levels) to the Dataset
    dataset.feat_ids = feat_ids
    dataset.feat_id_levels = feat_id_levels

    # close the database connection
    con.close()


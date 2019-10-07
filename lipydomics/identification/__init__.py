"""
    lipydomics/identification/__init__.py
    Dylan H. Ross
    2019/10/04

    description:
        Module for performing identification of individual lipid features
"""


from sqlite3 import connect


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

        The identifications are stored in `Dataset.ids` and the associated identification levels are stored in 
        `Dataset.id_levels`

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
    if level not in ['theo_mz', 'meas_mz_rt_ccs', 'any']:
        m = 'add_feature_ids: identification level "{}" not recognized'
        raise ValueError(m.format(level))

    # initialize connection to lipids.db

    ids, id_levels = [], []
    for mz, rt, ccs in dataset.labels:

        # try other identification strategies ...

        # fallback identification
        feat_id = 'UNK_{:09.4f}_{:05.2f}_{:06.2f}'.format(mz, rt, ccs)
        ids.append(feat_id)
        id_levels.append('')

    # apply the identifications (and identification levels) to the Dataset
    dataset.ids = ids
    dataset.id_levels = id_levels


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


def id_feature_theo_mz(cursor, mz, tol_mz, esi_mode):
    """
id_feat_theo_mz
    description:
        identifies a feature on the basis of theoretical m/z
    parameters:
        cursor () -- cursor for querying lipids.db
        mz (float) -- m/z to match
        tol_mz (float) -- tolerance for m/z search
        esi_mode (str) -- filter results by ionization mode: 'neg', 'pos', or None for unspecified
    returns:
        (str or list(str)), (str) -- putative identification(s) (or '' for no matches), identification level
"""
    return '', ''


def id_feature_meas_mz_rt_ccs(cursor, mz, rt, ccs, tol_mz, tol_rt, tol_ccs, esi_mode):
    """
id_feat_theo_mz
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
    return '', ''

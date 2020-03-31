"""
    lipydomics/identification/__init__.py
    Dylan H. Ross
    2019/10/04

    description:
        Module for performing identification of individual lipid features
"""


from sqlite3 import connect
import os


from lipydomics.identification.id_levels import (
    id_feat_any, id_feat_meas_mz_rt_ccs, id_feat_theo_mz_rt_ccs, id_feat_meas_mz_rt, id_feat_theo_mz_rt,
    id_feat_meas_mz_ccs, id_feat_theo_mz_ccs, id_feat_meas_mz, id_feat_theo_mz
)


def add_feature_ids(dataset, tol, level='any', norm='l2', db_version_tstamp=None, use_rt=True):
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
            'theo_mz_rt' -- match on theoretical m/z and rt
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

        The same tolerances are used for all levels of identification, and CCS tolerance is a percentage NOT an absolute
        number

        Potential identifications are automatically restricted on the basis of electrospray ionization mode (taken from
        Dataset.esi_mode) if it is set to 'pos' or 'neg'

        > If a retention time calibration is available, use calibrated retention time for feature identification <

        * if identifications have already been made, subsequent calls to this function will override previous results *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        tol (list(float, float, float)) -- tolerance for m/z, rt, and CCS, respectively
        [level (str)] -- specify the level of identification [optional, default='all']
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
        [db_version_tstamp (str or None)] -- use a specific time-stamped version of the lipids database instead of the
                                             default (most recent build) [optional, default=None]
        [use_rt (bool)] -- whether to use identification levels that involve retention time, ignored unless used with
                            the 'any' identification level [optional, default=True]
"""
    if level not in ['theo_mz', 'theo_mz_ccs', 'theo_mz_rt_ccs', 'meas_mz_ccs', 'meas_mz_rt_ccs', 'any']:
        m = 'add_feature_ids: identification level "{}" not recognized'
        raise ValueError(m.format(level))

    if db_version_tstamp:  # option to use a time-stamped version from the builds directory
        db_path = os.path.join(os.path.dirname(__file__), 'builds/lipids_{}.db'.format(db_version_tstamp))
    else:
        db_path = os.path.join(os.path.dirname(__file__), 'lipids.db')
    # make sure we can find the database
    if not os.path.isfile(db_path):
        m = 'add_feature_ids: unable to find lipid database ({})'.format(db_path)
        raise RuntimeError(m)

    # available identification functions
    id_funcs = {
        'any': id_feat_any,
        'meas_mz_rt_ccs': id_feat_meas_mz_rt_ccs,
        'theo_mz_rt_ccs': id_feat_theo_mz_rt_ccs,
        'meas_mz_rt': id_feat_meas_mz_rt,
        'theo_mz_rt': id_feat_theo_mz_rt,
        'meas_mz_ccs': id_feat_meas_mz_ccs,
        'theo_mz_ccs': id_feat_theo_mz_ccs,
        'meas_mz': id_feat_meas_mz,
        'theo_mz': id_feat_theo_mz
    } 

    # ESI mode from Dataset
    esi = dataset.esi_mode

    # initialize connection to lipids.db (stored within the lipydomics package)
    con = connect(db_path)
    cur = con.cursor()

    feat_ids, feat_id_levels, feat_id_scores = [], [], []
    for mz, rt, ccs in dataset.labels:

        # use calibrated retention time if a retention time calibration has been set up
        rt = dataset.rt_calibration.get_calibrated_rt(rt) if dataset.rt_calibration is not None else rt

        # convert the CCS tolerance into an absolute from the percentage
        tol2 = tol
        tol2[2] = tol2[2] / 100. * ccs

        # try to get identification(s)
        if level == 'any' and not use_rt:
            # use any identification level that does not include retention time
            feat_id, feat_id_level, feat_id_score = id_funcs['any'](cur, mz, rt, ccs, *tol2, esi,
                                                                    norm=norm, use_rt=False)
        else:
            feat_id, feat_id_level, feat_id_score = id_funcs[level](cur, mz, rt, ccs, *tol2, esi, norm=norm)

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


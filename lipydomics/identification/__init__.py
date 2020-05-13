"""
    lipydomics/identification/__init__.py
    Dylan H. Ross
    2019/10/04

    description:
        Module for performing identification of individual lipid features
"""


from sqlite3 import connect
import os
import warnings
import pickle


from lipydomics.identification.id_levels import (
    id_feat_any, id_feat_meas_mz_rt_ccs, id_feat_theo_mz_rt_ccs, id_feat_meas_mz_rt, id_feat_theo_mz_rt,
    id_feat_meas_mz_ccs, id_feat_theo_mz_ccs, id_feat_meas_mz, id_feat_theo_mz, id_feat_custom
)
from lipydomics.identification.encoder_params import (
    ccs_lipid_classes, ccs_ms_adducts, ccs_fa_mods, rt_lipid_classes, rt_fa_mods
)
from lipydomics.identification.train_lipid_ccs_pred import (
    prep_encoders as ccs_prep_encoders, featurize as ccs_featurize
)
from lipydomics.identification.mz_generation import get_lipid_mz
from lipydomics.identification.train_lipid_rt_pred import (
    prep_encoders as rt_prep_encoders, featurize as rt_featurize
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
            'meas_mz' -- simple matching based on measured m/z
            'theo_mz_ccs' -- match on theoretical m/z and CCS\
            'meas_mz_ccs' -- match on measured m/z and CCS
            'theo_mz_rt' -- match on theoretical m/z and rt
            'meas_mz_rt' -- match on measured m/z and rt
            'theo_mz_rt_ccs' -- match on theoretical m/z, rt, and CCS
            'meas_mz_rt_ccs' -- match on m/z, rt, and CCS
            (high)
            'any' -- start at the highest level, then work downward

        The `level` param can also be a list of specific levels, which will be attempted in the order provided.

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
        [level (str or list(str))] -- specify a single level of confidence for identifications, a list of confidence
                                        levels, or 'any' to use a tiered approach [optional, default='any']
        [norm (str)] -- specify l1 or l2 norm for computing scores [optional, default='l2']
        [db_version_tstamp (str or None)] -- use a specific time-stamped version of the lipids database instead of the
                                             default (most recent build) [optional, default=None]
        [use_rt (bool)] -- whether to use identification levels that involve retention time, ignored unless used with
                            the 'any' identification level [optional, default=True]
"""
    if level not in ['theo_mz', 'theo_mz_ccs', 'theo_mz_rt_ccs', 'meas_mz_ccs', 'meas_mz_rt_ccs', 'any',
                     'meas_mz', 'meas_mz_rt', 'theo_mz_rt'] and type(level) is not list:
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
        with warnings.catch_warnings():  # ignore a warning that happens with this division
            warnings.simplefilter("ignore")
            tol2[2] = tol2[2] / 100. * ccs

        # try to get identification(s)
        if type(level) is list:
            # use custom list of identification levels
            feat_id, feat_id_level, feat_id_score = id_feat_custom(level, cur, mz, rt, ccs, *tol2, esi, norm=norm)

        elif level == 'any' and not use_rt:
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


def predict_ccs(lipid_class, lipid_nc, lipid_nu, adduct, mz='generate', fa_mod=None, ignore_encoding_errors=False):
    """
predict_ccs
    description:
        Predicts a theoretical CCS of a lipid as defined by lipid class, fatty acid sum composition and MS adduct. If
        any of the parameters are not specifically encoded (i.e. not present in the training data) a ValueError is
        raised. This behavior can be overridden by the ignore_encoding_errors flag in order to get the prediction to be
        made regardless, although such predictions are subject to a high degree of error.
    parameters:
        lipid_class (str) -- lipid class
        lipid_nc (int) -- sum composition, number of fatty acid carbons
        lipid_nu (int) -- sum composition, number of fatty acid unsaturations
        adduct (str) -- MS adduct
        [mz (str or float)] -- m/z of MS adduct or 'generate' to generate an m/z value automatically using LipidMass
                                [optional, default='generate']
        [fa_mod (None or str)] -- fatty acid modifier (e.g. 'p', 'o') [optional, default=None]
        [ignore_encoding_errors (bool)] -- generate a prediction even if one or more of the input parameters are not
                                            encodable [optional, default=False]
    returns:
        (float) -- predicted CCS
"""
    # first check whether the lipid class, MS adduct and FA mod are encodable
    lc_ok = lipid_class in ccs_lipid_classes
    ad_ok = adduct in ccs_ms_adducts
    fm_ok = fa_mod is None or fa_mod in ccs_fa_mods
    all_ok = lc_ok and ad_ok and fm_ok
    if not (all_ok or ignore_encoding_errors):  # either all of the checks were good or we are ignoring errors
        m = ''
        if not lc_ok:
            m += 'lipid class "{}" not encodable '.format(lipid_class)
        if not ad_ok:
            m += 'MS adduct "{}" not encodable '.format(adduct)
        if not fm_ok:
            m += 'FA modifier "{}" not encodable '.format(fa_mod)
        raise ValueError('predict_ccs: {}'.format(m))

    # try to generate an m/z value if one wasn't provided
    if mz == 'generate':
        try:
            mz = get_lipid_mz(lipid_class, lipid_nc, lipid_nu, adduct, fa_mod=fa_mod)
        except ValueError as ve:
            m = 'predict_ccs: unable to generate m/z for lipid: "{}({}{}:{})_{}" ({})'
            m = m.format(lipid_class, '' if fa_mod is None else fa_mod, lipid_nc, lipid_nu, adduct, ve)
            raise ValueError(m)

    # prepare encoders
    c_encoder, f_encoder, a_encoder = ccs_prep_encoders()

    # load the predictive model and the scaler
    this_dir = os.path.dirname(__file__)
    model_path = os.path.join(this_dir, 'lipid_ccs_pred.pickle')
    scaler_path = os.path.join(this_dir, 'lipid_ccs_scale.pickle')
    with open(model_path, 'rb') as pf1, open(scaler_path, 'rb') as pf2:
        model = pickle.load(pf1)
        scaler = pickle.load(pf2)

    # featurize, scale, and predict CCS
    x = [ccs_featurize(lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, c_encoder, f_encoder, a_encoder)]
    return model.predict(scaler.transform(x))[0]


def predict_rt(lipid_class, lipid_nc, lipid_nu, fa_mod=None, ignore_encoding_errors=False):
    """
predict_ccs
    description:
        Predicts a theoretical HILIC retention time of a lipid as defined by lipid class and fatty acid sum composition.
        If any of the parameters are not specifically encoded (i.e. not present in the training data) a ValueError is
        raised. This behavior can be overridden by the ignore_encoding_errors flag in order to get the prediction to be
        made regardless, although such predictions are subject to a high degree of error.
    parameters:
        lipid_class (str) -- lipid class
        lipid_nc (int) -- sum composition, number of fatty acid carbons
        lipid_nu (int) -- sum composition, number of fatty acid unsaturations
        [fa_mod (None or str)] -- fatty acid modifier (e.g. 'p', 'o') [optional, default=None]
        [ignore_encoding_errors (bool)] -- generate a prediction even if one or more of the input parameters are not
                                            encodable [optional, default=False]
    returns:
        (float) -- predicted HILIC retention time
"""
    # first check whether the lipid class and FA mod are encodable
    lc_ok = lipid_class in ccs_lipid_classes
    fm_ok = fa_mod is None or fa_mod in ccs_fa_mods
    all_ok = lc_ok and fm_ok
    if not (all_ok or ignore_encoding_errors):  # either all of the checks were good or we are ignoring errors
        m = ''
        if not lc_ok:
            m += 'lipid class "{}" not encodable '.format(lipid_class)
        if not fm_ok:
            m += 'FA modifier "{}" not encodable '.format(fa_mod)
        raise ValueError('predict_rt: {}'.format(m))

    # prepare encoders
    c_encoder, f_encoder = rt_prep_encoders()

    # load the predictive model and the scaler
    this_dir = os.path.dirname(__file__)
    model_path = os.path.join(this_dir, 'lipid_rt_pred.pickle')
    scaler_path = os.path.join(this_dir, 'lipid_rt_scale.pickle')
    with open(model_path, 'rb') as pf1, open(scaler_path, 'rb') as pf2:
        model = pickle.load(pf1)
        scaler = pickle.load(pf2)

    # featurize, scale, and predict RT
    x = [rt_featurize(lipid_class, lipid_nc, lipid_nu, fa_mod, c_encoder, f_encoder)]
    return model.predict(scaler.transform(x))[0]

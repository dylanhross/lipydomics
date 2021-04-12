"""
    lipydomics/identification/train_lipid_ccs_pred.py
    Dylan H. Ross
    2019/10/08

    description:
        Trains a predictive model for generating predicted CCS values

        * requires scikit-learn v0.21.3 ! *
"""


from sqlite3 import connect
import os
import pickle
import numpy as np
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

from ..util import print_and_log
from .build_params import ccs_pred_ref_dsets
from .encoder_params import ccs_lipid_classes, ccs_fa_mods, ccs_ms_adducts


def prep_encoders():
    """
prep_encoders
    description:
        fits and returns encoders for lipid_class, fa_mod, and adduct
    returns:
        c_encoder, f_encoder, a_encoder (sklearn.preprocessing.OneHotEncoder) -- encoders for lipid_class, fa_mod, and
                                                                                    adduct, respectively
"""
    lipid_classes = [[_] for _ in ccs_lipid_classes]
    c_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(lipid_classes)
    fa_mods = [[_] for _ in ccs_fa_mods]
    f_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(fa_mods)
    adducts = [[_] for _ in ccs_ms_adducts]
    a_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(adducts)
    return c_encoder, f_encoder, a_encoder


def featurize(lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, c_encoder, f_encoder, a_encoder):
    """
featurize
    description:
        generates a numerical representation for a given lipid
    parameters:
        lipid_class (str) -- lipid class
        lipid_nc (int) -- sum composition: number of carbons
        lipid_nu (int) -- sum composition: number of unsaturations
        fa_mod (str) -- fatty acid modifiers
        adduct (str) -- MS adduct
        mz (float) -- m/z
        c_encoder, f_encoder, a_encoder (sklearn.preprocessing.OneHotEncoder) -- encoders for lipid_class, fa_mod, and
                                                                                    adduct, respectively
    returns:
        (np.array(float)) -- feature vector
"""
    lc_enc = c_encoder.transform([[lipid_class]])[0]
    fm_enc = f_encoder.transform([[fa_mod]])[0]
    ad_enc = a_encoder.transform([[adduct]])[0]
    lnc = np.array([float(lipid_nc)])
    lnu = np.array([float(lipid_nu)])
    m = np.array([mz])
    return np.concatenate([lc_enc, fm_enc, ad_enc, lnc, lnu, m])


def train_new_model(cursor, use_model, bl):
    """
train_new_model
    description:
        trains a predictive model
    parameters:
        cursor (sqlite3.cursor) -- cursor for querying lipids.db
        use_model (str) -- specify the ML model to use
        bl (file) -- build log
    returns:
        mdl, scaler -- trained predictive model and input scaler instances
"""
    # prepare encoders
    c_encoder, f_encoder, a_encoder = prep_encoders()

    # report the datasets used to train the predictive model
    print_and_log('using these datasets for CCS prediction:\n\t', bl, end='')
    print_and_log(' '.join(ccs_pred_ref_dsets), bl)

    # get the raw data and featurize (encode lipid_class, fa_mod, and adduct)
    qry = 'SELECT lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, ccs, src_tag FROM measured'
    X, y = [], []
    for lc, lnc, lnu, fam, add, m, c, st in cursor.execute(qry).fetchall():
        # only use the classes, fa_mods and adducts that are explicitly encoded
        lc_ok = lc in ccs_lipid_classes
        fam_ok = fam is None or fam in ccs_fa_mods
        add_ok = add in ccs_ms_adducts
        # only include data from the designated sources
        src_ok = st in ccs_pred_ref_dsets
        if lc_ok and fam_ok and add_ok and src_ok:
            X.append(featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder))
            y.append(float(c))
    X, y = np.array(X), np.array(y)
    print_and_log('X: {}'.format(X.shape), bl)
    print_and_log('y: {}'.format(y.shape), bl)

    # split into test/train sets, scale data (do not center)
    print_and_log('splitting data into training and test sets', bl)
    SSplit = ShuffleSplit(n_splits=1, test_size=0.2, random_state=1234)
    for train_index, test_index in SSplit.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
    print_and_log('X_train: {}'.format(X_train.shape), bl)
    print_and_log('y_train: {}'.format(y_train.shape), bl)
    print_and_log('X_test: {}'.format(X_test.shape), bl)
    print_and_log('y_test: {}'.format(y_test.shape), bl)
    print_and_log('scaling input data', bl)
    scaler = StandardScaler(with_mean=False)
    X_train_s = scaler.fit_transform(X_train)
    print_and_log('X_train_s: {}'.format(X_train_s.shape), bl)

    # train model (hyperparameters have already been tuned using gridsearch)
    print_and_log('training model', bl)
    if use_model == 'linear':
        model = LinearRegression(n_jobs=-1)
    elif use_model == 'svr':
        model = SVR(kernel='rbf', cache_size=1024, gamma='scale', C=4000)
    elif use_model == 'forest':
        model = RandomForestRegressor(n_jobs=-1, random_state=1234, bootstrap=True, max_depth=16, n_estimators=16)
    else:
        m = 'train_new_model: ML model "{}" not recognized'
        raise ValueError(m.format(use_model))
    model.fit(X_train_s, y_train)

    # performance on training set
    print_and_log('TRAINING SET PERFORMANCE', bl)
    y_train_pred = model.predict(X_train_s)
    y_train_abs_err = np.abs(y_train_pred - y_train)
    y_train_rel_err = 100. * (y_train_abs_err / y_train)
    print_and_log('mean absolute error: {:.2f} Å^2'.format(np.mean(y_train_abs_err)), bl)
    print_and_log('median absolute error: {:.2f} Å^2'.format(np.median(y_train_abs_err)), bl)
    print_and_log('median relative error: {:.2f} %'.format(np.median(y_train_rel_err)), bl)
    print_and_log('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_train, y_train_pred))), bl)

    # performance on test set
    print_and_log('TEST SET PERFORMANCE', bl)
    y_test_pred = model.predict(scaler.transform(X_test))
    y_test_abs_err = np.abs(y_test_pred - y_test)
    y_test_rel_err = 100. * (y_test_abs_err / y_test)
    print_and_log('mean absolute error: {:.2f} Å^2'.format(np.mean(y_test_abs_err)), bl)
    print_and_log('median absolute error: {:.2f} Å^2'.format(np.median(y_test_abs_err)), bl)
    print_and_log('median relative error: {:.2f} %'.format(np.median(y_test_rel_err)), bl)
    print_and_log('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_test, y_test_pred))), bl)

    # save the model and the scaler
    this_dir = os.path.dirname(__file__)
    model_path = os.path.join(this_dir, 'lipid_ccs_pred.pickle')
    scaler_path = os.path.join(this_dir, 'lipid_ccs_scale.pickle')
    with open(model_path, 'wb') as pf1, open(scaler_path, 'wb') as pf2:
        pickle.dump(model, pf1)
        pickle.dump(scaler, pf2)

    # return model and scaler
    return model, scaler


def dump_split_data_to_files(savedir):
    """
dump_split_data_to_files
    description:
        assembles training/test datasets just as would be done for actual model training then dumps those into
        separate .csv files: 'train.csv' and 'test.csv'
    parameters:
        savedir (str) -- directory to save the dumped files into
"""
    # connect to database
    db_path = os.path.join(os.path.dirname(__file__), 'lipids.db')
    con = connect(db_path)
    cur = con.cursor()

    # prepare encoders
    c_encoder, f_encoder, a_encoder = prep_encoders()

    qry = 'SELECT lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, ccs, src_tag FROM measured'
    X, y, l = [], [], []
    for lc, lnc, lnu, fam, add, m, c, st in cur.execute(qry).fetchall():
        # only use the classes, fa_mods and adducts that are explicitly encoded
        lc_ok = lc in ccs_lipid_classes
        fam_ok = fam is None or fam in ccs_fa_mods
        add_ok = add in ccs_ms_adducts
        # only include data from the designated sources
        src_ok = st in ccs_pred_ref_dsets
        if lc_ok and fam_ok and add_ok and src_ok:
            X.append(featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder))
            y.append(float(c))
            l.append('{}({}{}:{})_{}'.format(lc, fam if fam is not None else '', lnc, lnu, add))
    X, y, l = np.array(X), np.array(y), np.array(l)

    # split into test/train sets, scale data (do not center)
    SSplit = ShuffleSplit(n_splits=1, test_size=0.2, random_state=1234)
    for train_index, test_index in SSplit.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        l_train, l_test = l[train_index], l[test_index]

    # dump each array to file
    np.savetxt(os.path.join(savedir, 'X_train.csv'), X_train, delimiter=',')
    np.savetxt(os.path.join(savedir, 'y_train.csv'), y_train, delimiter=',')
    np.savetxt(os.path.join(savedir, 'l_train.csv'), l_train, delimiter=',', fmt='%s')
    np.savetxt(os.path.join(savedir, 'X_test.csv'), X_test, delimiter=',')
    np.savetxt(os.path.join(savedir, 'y_test.csv'), y_test, delimiter=',')
    np.savetxt(os.path.join(savedir, 'l_test.csv'), l_test, delimiter=',', fmt='%s')

    # close db connection
    con.close()


def main(tstamp):
    """ main build function """

    # connect to database
    db_path = os.path.join(os.path.dirname(__file__), 'lipids.db')
    con = connect(db_path)
    cur = con.cursor()

    # prepare encoders
    c_encoder, f_encoder, a_encoder = prep_encoders()

    build_log = os.path.join(os.path.dirname(__file__), 'builds/build_log_{}.txt'.format(tstamp))
    with open(build_log, 'a') as bl:

        # train the new model
        print_and_log('training new predictive CCS model (and input scaler) ...', bl)
        model, scaler = train_new_model(cur, 'svr', bl)
        print('... ok')
        print('... ok', file=bl)

        # add predicted CCS to the database
        print_and_log('\nadding predicted CCS to database ...', bl, end=' ')
        qry = 'SELECT t_id, lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz FROM predicted_mz'
        tid_to_ccs = {}
        for tid, lc, lnc, lnu, fam, add, m in cur.execute(qry).fetchall():
            if int(sum(c_encoder.transform([[lc]])[0])) != 0:  # make sure lipid class is encodable
                x = [featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder)]
                tid_to_ccs[int(tid)] = model.predict(scaler.transform(x))[0]
        qry = 'INSERT INTO predicted_ccs VALUES (?, ?)'
        for tid in tid_to_ccs:
            cur.execute(qry, (tid, tid_to_ccs[tid]))
        print_and_log('ok\n', bl)

    # commit changes to the database and close connection
    con.commit()
    con.close()

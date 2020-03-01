"""
    lipydomics/identification/train_lipid_ccs_pred.py
    Dylan H. Ross
    2019/10/08

    description:
        Trains a predictive model for generating theoretical CCS values

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


def prep_encoders():
    """
prep_encoders
    description:
        fits and returns encoders for lipid_class, fa_mod, and adduct
    returns:
        c_encoder, f_encoder, a_encoder (sklearn.preprocessing.OneHotEncoder) -- encoders for lipid_class, fa_mod, and
                                                                                    adduct, respectively
"""
    lipid_classes = [
        ['PE'],
        ['PC'],
        ['PG'],
        ['PS'],
        ['SM'],
        ['GlcCer'],
        ['Cer'],
        ['PI'],
        ['LPE'],
        ['DGDG'],
        ['PA'],
        ['TG'],
        ['CL'],
        ['LPC'],
        ['AcylPG'],
        ['LysylPG'],
        ['MGDG'],
        ['LPG'],
        ['AcylPE'],
        ['DG'],
        ['LPS'],
        ['GlcADG'],
        ['LPI'],
        ['AlaPG'],
        ['PIP']
    ]
    c_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(lipid_classes)
    fa_mods = [['o'], ['p']]
    f_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(fa_mods)
    adducts = [
        ['[M+H]+'],
        ['[M-H]-'],
        ['[M+Na]+'],
        ['[M+NH4]+'],
        ['[M+2Na-H]+'],
        ['[M+HCOO]-'],
        ['[M+Cl]-'],
        ['[M+K]+'],
        ['[M+H-H2O]+'],
        ['[M-2H]2-'],
        ['[M+CH3COO]-'],
        ['[M+2K]2+'],
        ['[M+H-2H2O]+'],
        ['[M+Na-2H2O]+'],
        ['[M+Na-H2O]+']
    ]
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

    # get the raw data and featurize (encode lipid_class, fa_mod, and adduct)
    qry = 'SELECT lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, ccs FROM measured'
    X, y = [], []
    for lc, lnc, lnu, fam, add, m, c in cursor.execute(qry).fetchall():
        X.append(featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder))
        y.append(float(c))
    X, y = np.array(X), np.array(y)
    print('X: ', X.shape)
    print('y: ', y.shape)
    print('X: ', X.shape, file=bl)
    print('y: ', y.shape, file=bl)

    # split into test/train sets, scale data (do not center)
    print('splitting data into training and test sets')
    SSplit = ShuffleSplit(n_splits=1, test_size=0.2, random_state=1234)
    for train_index, test_index in SSplit.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
    print('X_train: ', X_train.shape)
    print('y_train: ', y_train.shape)
    print('X_test: ', X_test.shape)
    print('y_test: ', y_test.shape)
    print('X_train: ', X_train.shape, file=bl)
    print('y_train: ', y_train.shape, file=bl)
    print('X_test: ', X_test.shape, file=bl)
    print('y_test: ', y_test.shape, file=bl)
    print('scaling input data')
    print('scaling input data', file=bl)
    scaler = StandardScaler(with_mean=False)
    X_train_s = scaler.fit_transform(X_train)
    print('X_train_s: ', X_train_s.shape)
    print('X_train_s: ', X_train_s.shape, file=bl)

    # train model (hyperparameters have already been tuned using gridsearch)
    print('training model')
    print('training model', file=bl)
    if use_model == 'linear':
        model = LinearRegression(n_jobs=-1)
    elif use_model == 'svr':
        model = SVR(kernel='rbf', cache_size=1024, gamma='scale', C=3000)
    elif use_model == 'forest':
        model = RandomForestRegressor(n_jobs=-1, random_state=1234, bootstrap=True, max_depth=16, n_estimators=8)
    else:
        m = 'train_new_model: ML model "{}" not recognized'
        raise ValueError(m.format(use_model))
    model.fit(X_train_s, y_train)

    # performance on training set
    print('TRAINING SET PERFORMANCE')
    print('TRAINING SET PERFORMANCE', file=bl)
    y_train_pred = model.predict(X_train_s)
    y_train_abs_err = np.abs(y_train_pred - y_train)
    print('mean absolute error: {:.2f} Å^2'.format(np.mean(y_train_abs_err)))
    print('median absolute error: {:.2f} Å^2'.format(np.median(y_train_abs_err)))
    print('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_train, y_train_pred))))
    print('mean absolute error: {:.2f} Å^2'.format(np.mean(y_train_abs_err)), file=bl)
    print('median absolute error: {:.2f} Å^2'.format(np.median(y_train_abs_err)), file=bl)
    print('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_train, y_train_pred))), file=bl)

    # performance on test set 
    print('TEST SET PERFORMANCE')
    print('TEST SET PERFORMANCE', file=bl)
    y_test_pred = model.predict(scaler.transform(X_test))
    y_test_abs_err = np.abs(y_test_pred - y_test)
    print('mean absolute error: {:.2f} Å^2'.format(np.mean(y_test_abs_err)))
    print('median absolute error: {:.2f} Å^2'.format(np.median(y_test_abs_err)))
    print('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_test, y_test_pred))))
    print('mean absolute error: {:.2f} Å^2'.format(np.mean(y_test_abs_err)), file=bl)
    print('median absolute error: {:.2f} Å^2'.format(np.median(y_test_abs_err)), file=bl)
    print('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_test, y_test_pred))), file=bl)

    # save the model and the scaler
    this_dir = os.path.dirname(__file__)
    model_path = os.path.join(this_dir, 'lipid_ccs_pred.pickle')
    scaler_path = os.path.join(this_dir, 'lipid_ccs_scale.pickle')
    with open(model_path, 'wb') as pf1, open(scaler_path, 'wb') as pf2:
        pickle.dump(model, pf1)
        pickle.dump(scaler, pf2)

    # return model and scaler
    return model, scaler


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
        print('training new predictive CCS model (and input scaler) ...')
        print('training new predictive CCS model (and input scaler) ...', file=bl)
        model, scaler = train_new_model(cur, 'svr', bl)
        print('... ok')
        print('... ok', file=bl)

        # add theoretical CCS to the database
        print('\nadding predicted CCS to database ...', end=' ')
        print('\nadding predicted CCS to database ...', end=' ', file=bl)
        qry = 'SELECT t_id, lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz FROM theoretical_mz'
        tid_to_ccs = {}
        for tid, lc, lnc, lnu, fam, add, m in cur.execute(qry).fetchall():
            if int(sum(c_encoder.transform([[lc]])[0])) != 0:  # make sure lipid class is encodable
                x = [featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder)]
                tid_to_ccs[int(tid)] = model.predict(scaler.transform(x))[0]
        qry = 'INSERT INTO theoretical_ccs VALUES (?, ?)'
        for tid in tid_to_ccs:
            cur.execute(qry, (tid, tid_to_ccs[tid]))
        print('ok\n')
        print('ok\n', file=bl)

        # commit changes to the database and close connection
        con.commit()
        con.close()

#!/usr/local/Cellar/python3/3.7.3/bin/python3
"""
    lipydomics/identification/train_lipid_ccs_pred.py
    Dylan H. Ross
    2019/10/08

    description:
        Trains a predictive model for generating theoretical CCS values

        * requires scikit-learn v0.21.3 ! *
"""


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
        ['AcylPG'],
        ['LPG'],
        ['LysylPG'],
        ['LPC'],
        ['CL'],
        ['PA'],
        ['TG'],
        ['SM'],
        ['GlcCer'],
        ['Cer'],
        ['DGDG'],
        ['LPE'],
        ['PI'],
        ['PS'],
        ['PC'],
        ['PG'],
        ['PE']
    ]
    c_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(lipid_classes)
    fa_mods = [['o'], ['p']]
    f_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(fa_mods)
    adducts = [
        ['[M+H-H2O]+'],
        ['[M+Cl]-'],
        ['[M+2K]2+'],
        ['[M+CH3COO]-'],
        ['[M-2H]2-'],
        ['[M+HCOO]-'],
        ['[M+K]+'],
        ['[M+NH4]+'],
        ['[M+H]+'],
        ['[M+Na]+'],
        ['[M-H]-'],
    ]
    a_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(adducts)
    return c_encoder, f_encoder, a_encoder


def featurize(lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, c_encoder, f_encoder, a_encoder):
    """
featurize
    description:

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


def train_new_model(cursor, use_model):
    """
train_new_model
    description:
        trains a predictive model
    parameters:
        cursor (sqlite3.cursor) -- cursor for querying lipids.db
        use_model (str) -- specify the ML model to use
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
    print('scaling input data')
    scaler = StandardScaler(with_mean=False)
    X_train_s = scaler.fit_transform(X_train)
    print('X_train_s: ', X_train_s.shape)

    # train model (hyperparameters have already been tuned using gridsearch)
    print('training model')
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
    y_train_pred = model.predict(X_train_s)
    y_train_abs_err = np.abs(y_train_pred - y_train)
    print('mean absolute error: {:.2f} Å^2'.format(np.mean(y_train_abs_err)))
    print('median absolute error: {:.2f} Å^2'.format(np.median(y_train_abs_err)))
    print('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_train, y_train_pred))))

    # performance on test set 
    print('TEST SET PERFORMANCE')
    y_test_pred = model.predict(scaler.transform(X_test))
    y_test_abs_err = np.abs(y_test_pred - y_test)
    print('mean absolute error: {:.2f} Å^2'.format(np.mean(y_test_abs_err)))
    print('median absolute error: {:.2f} Å^2'.format(np.median(y_test_abs_err)))
    print('RMSE: {:.2f} Å^2'.format(np.sqrt(mean_squared_error(y_test, y_test_pred))))

    # save the model and the scaler
    with open('lipid_ccs_pred.pickle', 'wb') as pf1, open('lipid_ccs_scale.pickle', 'wb') as pf2:
        pickle.dump(model, pf1)
        pickle.dump(scaler, pf2)

    # return model and scaler
    return model, scaler


if __name__ == '__main__':

    from sqlite3 import connect
    import os

    # connect to the database
    con = connect('lipids.db')
    cur = con.cursor()

    # prepare encoders
    c_encoder, f_encoder, a_encoder = prep_encoders() 
    
    # load the trained model if available, or train a new one
    model_path = 'lipid_ccs_pred.pickle'
    scaler_path = 'lipid_ccs_scale.pickle' 
    if not os.path.isfile(model_path) or not os.path.isfile(scaler_path):
        print('training new predictive CCS model (and input scaler) ...',)
        model, scaler = train_new_model(cur, 'svr')
        print('... ok')
    else:
        print('loading pre-trained predictive CCS model (and input scaler) ...', end=' ')
        with open(model_path, 'rb') as pf1, open(scaler_path, 'rb') as pf2:
            model = pickle.load(pf1)
            scaler = pickle.load(pf2)
        print('ok')

    # add theoretical CCS to the database
    print('adding predicted CCS to database ...', end=' ')
    qry = 'SELECT t_id, lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz FROM theoretical_mz'
    tid_to_ccs = {}
    for tid, lc, lnc, lnu, fam, add, m in cur.execute(qry).fetchall():

        x = [featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder)]
        tid_to_ccs[int(tid)] = model.predict(scaler.transform(x))[0]
    qry = 'INSERT INTO theoretical_ccs VALUES (?, ?)'
    for tid in tid_to_ccs:
        cur.execute(qry, (tid, tid_to_ccs[tid]))
    print('ok')

    # commit changes to the database and close connection
    con.commit()
    con.close()

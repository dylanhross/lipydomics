#!/usr/local/Cellar/python3/3.7.3/bin/python3
"""
    lipydomics/identification/train_lipid_rt_pred.py
    Dylan H. Ross
    2019/12/03

    description:
        Trains a predictive model for generating theoretical RT values

        * requires scikit-learn v0.21.3 ! *
"""

import pickle
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# reuse some utility functions from the CCS predictor training script
from train_lipid_ccs_pred import prep_encoders, featurize


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
    qry = 'SELECT lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz, rt FROM measured WHERE rt IS NOT NULL'
    X, y = [], []
    for lc, lnc, lnu, fam, add, m, c in cursor.execute(qry).fetchall():
        X.append(featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder))
        y.append(float(c))
    X, y = np.array(X), np.array(y)
    print('X: ', X.shape)
    print('y: ', y.shape)

    # split into test/train sets, scale data (do not center)
    print('splitting data into training and test sets')
    SSplit = ShuffleSplit(n_splits=1, test_size=0.2, random_state=1235)
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
        model = SVR(kernel='rbf', cache_size=1024, gamma='scale', C=300)
    elif use_model == 'forest':
        model = RandomForestRegressor(n_jobs=-1, random_state=1234, bootstrap=True, max_depth=8, n_estimators=4)
    else:
        m = 'train_new_model: ML model "{}" not recognized'
        raise ValueError(m.format(use_model))
    model.fit(X_train_s, y_train)

    # performance on training set
    print('TRAINING SET PERFORMANCE')
    y_train_pred = model.predict(X_train_s)
    y_train_abs_err = np.abs(y_train_pred - y_train)
    print('mean absolute error: {:.2f} min'.format(np.mean(y_train_abs_err)))
    print('median absolute error: {:.2f} min'.format(np.median(y_train_abs_err)))
    print('RMSE: {:.2f} min'.format(np.sqrt(mean_squared_error(y_train, y_train_pred))))

    # performance on test set
    print('TEST SET PERFORMANCE')
    y_test_pred = model.predict(scaler.transform(X_test))
    y_test_abs_err = np.abs(y_test_pred - y_test)
    print('mean absolute error: {:.2f} min'.format(np.mean(y_test_abs_err)))
    print('median absolute error: {:.2f} min'.format(np.median(y_test_abs_err)))
    print('RMSE: {:.2f} min'.format(np.sqrt(mean_squared_error(y_test, y_test_pred))))

    # save the model and the scaler
    with open('lipid_rt_pred.pickle', 'wb') as pf1, open('lipid_rt_scale.pickle', 'wb') as pf2:
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
    model_path = 'lipid_rt_pred.pickle'
    scaler_path = 'lipid_rt_scale.pickle'
    if not os.path.isfile(model_path) or not os.path.isfile(scaler_path):
        print('training new predictive RT model (and input scaler) ...', )
        model, scaler = train_new_model(cur, 'linear')
        print('... ok')
    else:
        print('loading pre-trained predictive RT model (and input scaler) ...', end=' ')
        with open(model_path, 'rb') as pf1, open(scaler_path, 'rb') as pf2:
            model = pickle.load(pf1)
            scaler = pickle.load(pf2)
        print('ok')

    # add theoretical CCS to the database
    print('\nadding predicted RT to database ...', end=' ')
    qry = 'SELECT t_id, lipid_class, lipid_nc, lipid_nu, fa_mod, adduct, mz FROM theoretical_mz'
    tid_to_rt = {}
    for tid, lc, lnc, lnu, fam, add, m in cur.execute(qry).fetchall():
        x = [featurize(lc, lnc, lnu, fam, add, float(m), c_encoder, f_encoder, a_encoder)]
        tid_to_rt[int(tid)] = model.predict(scaler.transform(x))[0]
    qry = 'INSERT INTO theoretical_rt VALUES (?, ?)'
    for tid in tid_to_rt:
        cur.execute(qry, (tid, tid_to_rt[tid]))
    print('ok\n')

    # commit changes to the database and close connection
    con.commit()
    con.close()

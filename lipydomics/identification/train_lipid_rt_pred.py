"""
    lipydomics/identification/train_lipid_rt_pred.py
    Dylan H. Ross
    2019/12/03

    description:
        Trains a predictive model for generating theoretical RT values

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
        c_encoder, f_encoder (sklearn.preprocessing.OneHotEncoder) -- encoders for lipid_class and fa_mod
"""
    lipid_classes = [
        ['PG'],
        ['PE'],
        ['DGDG'],
        ['LPE'],
        ['PC'],
        ['CL'],
        ['PI'],
        ['PA'],
        ['Cer'],
        ['AcylPG'],
        ['LysylPG'],
        ['GlcCer'],
        ['MGDG'],
        ['LPG'],
        ['AcylPE'],
        ['SM'],
        ['LPC'],
        ['PS'],
        ['DG'],
        ['GlcADG'],
        ['AlaPG'],
        ['PIP']
    ]
    c_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(lipid_classes)
    fa_mods = [['p']]
    f_encoder = OneHotEncoder(sparse=False, handle_unknown='ignore').fit(fa_mods)
    return c_encoder, f_encoder


def featurize(lipid_class, lipid_nc, lipid_nu, fa_mod, c_encoder, f_encoder):
    """
featurize
    description:
        generates a numerical representation for a given lipid
    parameters:
        lipid_class (str) -- lipid class
        lipid_nc (int) -- sum composition: number of carbons
        lipid_nu (int) -- sum composition: number of unsaturations
        fa_mod (str) -- fatty acid modifiers
        mz (float) -- m/z
        c_encoder, f_encoder (sklearn.preprocessing.OneHotEncoder) -- encoders for lipid_class and fa_mod
    returns:
        (np.array(float)) -- feature vector
"""
    lc_enc = c_encoder.transform([[lipid_class]])[0]
    fm_enc = f_encoder.transform([[fa_mod]])[0]
    lnc = np.array([float(lipid_nc)])
    lnu = np.array([float(lipid_nu)])
    return np.concatenate([lc_enc, fm_enc, lnc, lnu])


def train_new_model(cursor):
    """
train_new_model
    description:
        trains a predictive model
    parameters:
        cursor (sqlite3.cursor) -- cursor for querying lipids.db
    returns:
        mdl, scaler -- trained predictive model and input scaler instances
"""
    # prepare encoders
    c_encoder, f_encoder = prep_encoders()

    # get the raw data and featurize (encode lipid_class, fa_mod, and adduct)
    qry = 'SELECT lipid_class, lipid_nc, lipid_nu, fa_mod, rt FROM measured WHERE rt IS NOT NULL'
    X, y = [], []
    for lc, lnc, lnu, fam, c in cursor.execute(qry).fetchall():
        X.append(featurize(lc, lnc, lnu, fam, c_encoder, f_encoder))
        y.append(float(c))
    X, y = np.array(X), np.array(y)
    print('X: ', X.shape)
    print('y: ', y.shape)

    # split into test/train sets, scale data (do not center)
    print('splitting data into training and test sets')
    SSplit = ShuffleSplit(n_splits=1, test_size=0.2, random_state=1236)
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

    # train model
    print('training model')
    model = LinearRegression(n_jobs=-1)
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
    this_dir = os.path.dirname(__file__)
    model_path = os.path.join(this_dir, 'lipid_rt_pred.pickle')
    scaler_path = os.path.join(this_dir, 'lipid_rt_scale.pickle')
    with open(model_path, 'wb') as pf1, open(scaler_path, 'wb') as pf2:
        pickle.dump(model, pf1)
        pickle.dump(scaler, pf2)

    # return model and scaler
    return model, scaler


def main():
    """ main build function """

    # connect to database
    db_path = os.path.join(os.path.dirname(__file__), 'lipids.db')
    con = connect(db_path)
    cur = con.cursor()

    # prepare encoders
    c_encoder, f_encoder = prep_encoders()

    # train a new model
    print('training new predictive RT model (and input scaler) ...', )
    model, scaler = train_new_model(cur)
    print('... ok')

    # add theoretical CCS to the database
    print('\nadding predicted RT to database ...', end=' ')
    qry = 'SELECT t_id, lipid_class, lipid_nc, lipid_nu, fa_mod FROM theoretical_mz'
    tid_to_rt = {}
    for tid, lc, lnc, lnu, fam in cur.execute(qry).fetchall():
        if int(sum(c_encoder.transform([[lc]])[0])) != 0: # make sure lipid class is encodable
            x = [featurize(lc, lnc, lnu, fam, c_encoder, f_encoder)]
            tid_to_rt[int(tid)] = model.predict(scaler.transform(x))[0]
    qry = 'INSERT INTO theoretical_rt VALUES (?, ?)'
    for tid in tid_to_rt:
        cur.execute(qry, (tid, tid_to_rt[tid]))
    print('ok\n')

    # commit changes to the database and close connection
    con.commit()
    con.close()

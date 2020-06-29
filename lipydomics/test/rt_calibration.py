"""
    lipydomics/test/rt_calibration.py
    Dylan H. Ross
    2020/01/09

    description:
        Define tests associated with the lipydomics/identification/rt_calibration module
"""


import os
from csv import reader
from numpy import array, mean, sqrt

from lipydomics.test import run_tests
from lipydomics.identification.rt_calibration import get_ref_rt, RTCalibration
from lipydomics.util import parse_lipid


def get_ref_rt_lipids1():
    """
get_ref_rt_real1
    description:
        Tests the get_ref_rt helper function using the lipids from lipids_1.txt using both strict and non-strict modes

        Test fails if there are any errors
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    with open(os.path.join(os.path.dirname(__file__), 'lipids_1.txt'), 'r') as f:
        for line in f.readlines():
            parsed = parse_lipid(line.strip())
            fam = parsed['fa_mod'] if 'fa_mod' in parsed else None
            rt_strict = get_ref_rt(parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat'], fa_mod=fam)
            rt = get_ref_rt(parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat'], fa_mod=fam, strict=False)

    return True


def rtcal_init_mismatch_len():
    """
rtcal_init_mismatch_len
    description:
        Tests initializing an RTCalibration object with lists of mismatched length which triggers a ValueError

        Test fails if there are no errors or other errors besides the ValueError
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    try:
        rtc = RTCalibration(['A', 'B', 'C'], [1.23, 2.34, 3.45], [1.1, 2.2, 4.4, 8.8])
    except ValueError:
        return True
    return False


def rtcal_calibrate_rtc1_c12():
    """
rtcal_calibrate_rtc1_c12
    description:
        Uses the data from rtc_1.csv column 1 (reference rt) and column 2 (measured rt) to compute a RT calibration.
        Checks the root mean squared error of predicted vs. reference RT
        Uses 2 lipids for calibration

        Test fails if there are any errors or if RMSE of predicted RT is > 0.3 min
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    n = ['PG(36:0)', 'LysylPG(34:0)']
    rtm = [0.707, 3.928]
    rtr = [1.293, 6.598]
    rtc = RTCalibration(n, rtm, rtr)
    x, r = [], []
    with open(os.path.join(os.path.dirname(__file__), 'rtc_1.csv'), 'r') as f:
        rdr = reader(f)
        for l, rr, rm, _ in rdr:
            x.append(float(rm))
            r.append(float(rr))
    y = [rtc.get_calibrated_rt(_) for _ in x]
    y = array(y)
    r = array(r)
    rmse = sqrt(mean((y - r)**2.))
    if rmse > 0.3:
        return False
    return True


def rtcal_calibrate_rtc1_c13():
    """
rtcal_calibrate_rtc1_c13
    description:
        Uses the data from rtc_1.csv column 1 (reference rt) and column 3 (measured rt) to compute a RT calibration.
        Checks the root mean squared error of predicted vs. reference RT
        Uses 2 lipids for calibration

        Test fails if there are any errors or if RMSE of predicted RT is > 0.3 min
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    n = ['PG(36:0)', 'LysylPG(34:0)']
    rtm = [0.464, 2.533]
    rtr = [1.293, 6.598]
    rtc = RTCalibration(n, rtm, rtr)
    x, r = [], []
    with open(os.path.join(os.path.dirname(__file__), 'rtc_1.csv'), 'r') as f:
        rdr = reader(f)
        for l, rr, _, rm in rdr:
            x.append(float(rm))
            r.append(float(rr))
    y = [rtc.get_calibrated_rt(_) for _ in x]
    y = array(y)
    r = array(r)
    rmse = sqrt(mean((y - r)**2.))
    if rmse > 0.3:
        return False
    return True


def rtcal_calibrate_rtc2_c12():
    """
rtcal_calibrate_rtc1_c12
    description:
        Uses the data from rtc_2.csv column 1 (reference rt) and column 2 (measured rt) to compute a RT calibration.
        Checks the root mean squared error of predicted vs. reference RT
        Uses 5 lipids for calibration

        Test fails if there are any errors or if RMSE of predicted RT is > 0.2 min
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    n = ['FFA(18:0)', 'DGDG(36:0)', 'PG(36:0)', 'CL(60:0)', 'LysylPG(34:0)']
    rtm = [0.673, 0.929, 1.549, 2.843, 3.996]
    rtr = [0.673, 0.895, 1.549, 4.393, 7.252]
    rtc = RTCalibration(n, rtm, rtr)
    x, r = [], []
    with open(os.path.join(os.path.dirname(__file__), 'rtc_2.csv'), 'r') as f:
        rdr = reader(f)
        for l, rr, rm, _ in rdr:
            x.append(float(rm))
            r.append(float(rr))
    y = [rtc.get_calibrated_rt(_) for _ in x]
    y = array(y)
    r = array(r)
    rmse = sqrt(mean((y - r)**2.))
    if rmse > 0.2:
        return False
    return True


def rtcal_calibrate_rtc2_c13():
    """
rtcal_calibrate_rtc1_c13
    description:
        Uses the data from rtc_2.csv column 1 (reference rt) and column 3 (measured rt) to compute a RT calibration.
        Checks the root mean squared error of predicted vs. reference RT
        Uses 5 lipids for calibration

        Test fails if there are any errors or if RMSE of predicted RT is > 0.2 min
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    n = ['FFA(18:0)', 'DGDG(36:0)', 'PG(36:0)', 'CL(60:0)', 'LysylPG(34:0)']
    rtm = [0.673, 0.929, 1.516, 6.665, 12.246]
    rtr = [0.673, 0.895, 1.549, 4.393, 7.252]
    rtc = RTCalibration(n, rtm, rtr)
    x, r = [], []
    with open(os.path.join(os.path.dirname(__file__), 'rtc_2.csv'), 'r') as f:
        rdr = reader(f)
        for l, rr, _, rm in rdr:
            x.append(float(rm))
            r.append(float(rr))
    y = [rtc.get_calibrated_rt(_) for _ in x]
    y = array(y)
    r = array(r)
    rmse = sqrt(mean((y - r)**2.))
    if rmse > 0.2:
        return False
    return True


# references to al of the test functions to be run, and order to run them in
all_tests = [
    get_ref_rt_lipids1,
    rtcal_init_mismatch_len,
    rtcal_calibrate_rtc1_c12,
    rtcal_calibrate_rtc1_c13,
    rtcal_calibrate_rtc2_c12,
    rtcal_calibrate_rtc2_c13
]
if __name__ == '__main__':
    run_tests(all_tests)

"""
    lipydomics/identification/rt_calibration.py
    Dylan H. Ross
    2020/01/09

    description:
        Utilities for calibration of HILIC retention times that allow for comparison of retention times from different
        gradients.
"""


from sqlite3 import connect
import os
from numpy import mean


def get_ref_rt(lc, nc, nu, fa_mod=None, strict=True):
    """
get_ref_rt
    description:
        A helper function that searches the reference lipids database to find a reference HILIC retention time for
        a lipid using its class, number of carbons, and number of unsaturations. The number of unsaturations may
        optionally be ignored by setting the strict kwarg to False. If multiple reference entries match the search
        criteria, then the average retention time is returned. If no entries match the search criteria, then None is
        returned.
    parameters:
        lc (str) -- lipid class
        nc (int) -- sum composition number of carbons
        nu (int) -- sum composition number of unsaturations
        [fa_mod (str)] -- fatty acid modifier or None [optional, default=None]
        [strict (bool)] -- if False, allow a match only using lipid class and number of carbons [optional, default=True]
    returns:
        (float or None) -- retention time or None if unable to find a match
"""
    # initialize connection to lipids.db (stored within the lipydomics package)
    con = connect(os.path.join(os.path.dirname(__file__), 'lipids.db'))
    cur = con.cursor()
    nu_option = " AND lipid_nu={} ".format(nu) if strict else " "
    fa_mod_option = ' AND fa_mod="{}"'.format(fa_mod) if fa_mod else " AND fa_mod IS NULL "
    qry = "SELECT rt FROM measured WHERE lipid_class=? AND lipid_nc=? {}" \
          "{} AND rt IS NOT NULL".format(nu_option, fa_mod_option)
    qdata = (lc, nc)
    rts = [float(_[0]) for _ in cur.execute(qry, qdata).fetchall()]
    if rts == []:
        return None
    return mean(rts)


class RTCalibration:
    """
RTCalibration
    description:
        An object for performing HILIC retention time calibration
"""

    def __init__(self, lipids, meas_rt, ref_rt):
        """
RTCalibration.__init__
    description:
        Stores lists of lipid calibrants and their reference/measured retention times,
        lists are sorted together by measured retention time
    parameters:
        lipids (list(str)) -- lipid calibrants
        meas_rt (list(float)) -- measured retention times
        ref_rt (list(float)) -- reference retention times
"""
        # all the lists need to be the same length
        if len(lipids) != len(ref_rt) or len(lipids) != len(meas_rt):
            m = 'RTCalibration: __init__: lipids, ref_rt, and meas_rt must all be the same length ({}, {}, {})'
            raise ValueError(m.format(len(lipids), len(ref_rt), len(meas_rt)))
        # sort the calibrants by measured rt and store them
        self.meas_rt, self.ref_rt, self.lipids = (list(t) for t in zip(*sorted(zip(meas_rt, ref_rt, lipids))))
        self.n_calibrants = len(lipids)

    def __repr__(self):
        """
RTCalibration.__repr__
    description:
        generates a string representation of this RTCalibration instance
    returns:
        (str) -- string representation of this instance
"""
        s = 'RTCalibration(\n\tlipid, measured rt, reference rt\n'
        for lipid, ref, meas in zip(self.lipids, self.ref_rt, self.meas_rt):
            s += '\t{}, {:.2f} min, {:.2f} min\n'.format(lipid, meas, ref)
        s += ')'
        return s

    @staticmethod
    def linear(x1, y1, x2, y2, x):
        """
RTCalibration.linear
    description:
        Returns a linear interpolation (or extrapolation) of an x value using a line defined by two points
    parameters:
        x1 (float) -- point 1 x value
        y1 (float) -- point 1 y value
        x2 (float) -- point 2 x value
        y2 (float) -- point 2 y value
        x (float) -- input x value
    returns:
        (float) -- output y value
"""
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
        return m * x + b

    def get_calibrated_rt(self, rt_uncal):
        """
RTCalibration.get_calibrated_rt
    description:
        returns a calibrated retention time using the calibrant data
    parameters:
        rt_uncal (float) -- uncalibrated retention time
    returns:
         (float) -- calibrated retention time
"""
        # if there is only one calibrant, just compute a ratio and use that
        if self.n_calibrants == 1:
            return rt_uncal * (self.ref_rt[0] / self.meas_rt[0])

        # if there are two calibrants, use the line between them
        elif self.n_calibrants == 2:
            return self.linear(self.meas_rt[0], self.ref_rt[0], self.meas_rt[1], self.ref_rt[1], rt_uncal)

        # if there are more than two calibrants, use linear interpolation between each pair
        else:
            if rt_uncal <= self.meas_rt[1]:
                i0, i1 = 0, 1
            elif rt_uncal > self.meas_rt[-2]:
                i0, i1 = -2, -1
            else:
                i0, i1 = 1, 2
                while rt_uncal > self.meas_rt[i1]:
                    i0 += 1
                    i1 += 1
            return self.linear(self.meas_rt[i0], self.ref_rt[i0], self.meas_rt[i1], self.ref_rt[i1], rt_uncal)

        return None


def add_rt_calibration(dataset, lipids, meas_rt, ref_rt):
    """
add_feature_ids
    description:
        Creates an RTCalibration instance and stores it in a Dataset object

        The RTCalibration instance is stored in `Dataset.rt_calibration`

        * subsequent calls to this function will override previous results *
    parameters:
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        lipids (list(str)) -- lipid calibrants
        meas_rt (list(float)) -- measured retention times
        ref_rt (list(float)) -- reference retention times
"""
    dataset.rt_calibration = RTCalibration(lipids, meas_rt, ref_rt)

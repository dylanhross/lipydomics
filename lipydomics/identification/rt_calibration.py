"""
    lipydomics/identification/rt_calibration.py
    Dylan H. Ross
    2020/01/09

    description:
        Utilities for calibration of HILIC retention times that allow for comparison of retention times from different
        gradients.
"""


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

    def __init__(self):
        """
RTCalibration.__init__
    description:

    parameters:

"""
        pass

    def get_calibrated_rt(self, rt_uncal):
        """
RTCalibration
    description:
        returns a calibrated retention time using the fitted calibration function
    parameters:
        rt_uncal (float) -- uncalibrated retention time
    returns:
         (float) -- calibrated retention time
"""
        return None

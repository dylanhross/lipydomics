"""
    lipydomics/util.py
    Dylan H. Ross
    2020/03/17

    description:
        Miscellaneous utility functions that do not quite fit in with other modules
"""

import re
from numpy import sqrt, nan, inf, mean, abs, array
from datetime import datetime


def abbreviate_sheet(sheet_name):
    """
abbreviate_sheet
    description:
        For a pandas Dataframe to be exported into an Excel spreadsheet, the sheet names must be no longer than 31
        characters. In order to export the sheets containing computed statistics, the statistic labels often need to be
        shortened because they contain information on the statistic and the groups used to calculate it. This function
        attempts to abbreviate such labels to 31 characters. If a suitable abbreviation cannot be made, then the
        original name is simply truncated down to 31 characters. All group names are truncated to their first 3
        characters as well.
    parameters:
        sheet_name (str) -- sheet name
    returns:
        (str) -- abbreviated sheet name
"""
    # just in case this is called in error, always return the sheet name if it is already <32 characters
    if len(sheet_name) < 32:
        return sheet_name
    # abbreviations
    stat_abbrev = {'PCA3': 'PC3', 'PLS-DA': 'PLS', 'ANOVA': 'ANV', '2-group-corr': '2GC', 'LOG2FC': 'LFC'}
    load_or_proj_abbrev = {'': '', 'loadings': 'L', 'projections': 'P'}
    nrm_abbrev = {'raw': 'R', 'normed': 'N'}
    # parse the stat label for relevant info
    sheet_name_split = sheet_name.split('_')
    stat = sheet_name_split[0]
    if stat not in stat_abbrev:
        # fall back to simple truncation
        return sheet_name[:31]
    nrm = sheet_name_split[-1]
    load_or_proj = '' if sheet_name_split[-2] not in ['loadings', 'projections'] else sheet_name_split[-2]
    groups = '-'.join([name[:3] if len(name) > 3 else name for name in sheet_name_split[1].split('-')])
    # construct the abbreviated name
    abbrev = stat_abbrev[stat] + load_or_proj_abbrev[load_or_proj] + nrm_abbrev[nrm] + '_' + groups
    if len(abbrev) > 31:
        # fall back to simple truncation
        return sheet_name[:31]
    return abbrev


def filter_d(mzs, rts, ccss, data):
    """
filter_d
    description:
        a helper function for filtering data
        given M/Z, RT, CCS ranges and a DataFrame containing data,
        find and returns all data within that range.
        * CCS tolerance is absolute in this case, NOT a percentage *
    parameters:
        mzs (list(flaot)) -- mz [0] and tolerance [1]
        rts (list(float)) -- rt [0] and tolerance [1]
        ccss (list(float)) -- ccs [0] and tolerance [1]
        data (pandas.DataFrame) -- DataFrame representation of the data
"""
    # removed all of the casts to float, that should happen at input time not here
    filtered = data[(data[0] < mzs[0] + mzs[1]) & (data[0] > mzs[0] - mzs[1]) &
                    (data[1] < rts[0] + rts[1]) & (data[1] > rts[0] - rts[1]) &
                    (data[2] < ccss[0] + ccss[1]) & (data[2] > ccss[0] - ccss[1])]
    return filtered


def get_score(tol_mz, tol_rt, tol_ccs, mz_q=None, rt_q=None, ccs_q=None, mz_x=None, rt_x=None, ccs_x=None, norm='l2'):
    """
get_score
    description:
        computes a score reflecting the quality of an identification, using mz rt and ccs or any combination

        The score is determined by the residuals between the query values (q) and a potential match (x), normalized by
        their respective tolerances. If only a single pair of values (q and x) is provided, the score is simply the
        inverse of the normalized residual, otherwise, it is the inverse of the l1 or l2 norm of the normalized
        residuals vector. The norm kwarg controls whether the l1 or l2 norm is used in computing the score
    parameters:
        cursor (sqlite3.Cursor) -- cursor for querying lipids.db
        tol_mz (float) -- tolerance for m/z
        tol_rt (float) -- tolerance for retention time
        tol_ccs (float) -- tolerance for CCS
        [mz_q (None or float)] -- if specified, the query m/z [optional, default=None]
        [rt_q (None or float)] -- if specified, the query retention time [optional, default=None]
        [ccs_q (None or float)] -- if specified, the query CCS [optional, default=None]
        [mz_q (None or float)] -- if specified, the m/z of a potential match [optional, default=None]
        [rt_q (None or float)] -- if specified, the retention time of a potential match [optional, default=None]
        [ccs_q (None or float)] -- if specified, the CCS of a potential match [optional, default=None]
        [norm (str)] -- specify l1 or l2 norm [optional, default='l2']
    returns:
        (float) -- score (higher = more confidence in ID)
"""
    q = [mz_q, rt_q, ccs_q]
    x = [mz_x, rt_x, ccs_x]
    tol = [tol_mz, tol_rt, tol_ccs]
    # compute the normalized residuals
    rn = []
    for q_, x_, tol_ in zip(q, x, tol):
        if q_ and x_:
            rn.append((x_ - q_) / tol_)
    # if there was only a single value, just return the inverse of the absolute value
    if not rn:
        m = 'get_score: unable to compute residuals'
        raise RuntimeError(m)
    elif len(rn) == 1:
        return 1. / max(abs(rn[0]), 0.000001)  # prevent zero-division just in case ...
    else:
        if norm == 'l1':
            return 1. / max(sum([abs(_) for _ in rn]), 0.000001)  # prevent zero-division just in case ...
        elif norm == 'l2':
            return 1. / max(sqrt(sum([_**2. for _ in rn])), 0.000001)  # prevent zero-division just in case ...
        else:
            m = 'get_score: norm method "{}" not recognized'
            raise ValueError(m.format(norm))


def parse_lipid(name):
    """
parse_lipid
    description:
        parses a lipid name into lipid class and fatty acid composition, returning a
        dictionary with the information. Handles total fatty acid composition, as well
        as individual composition, examples:
            PC(38:3)        --> class: PC, n_carbon: 38, n_unsat: 3
            PC(18:1/20:2)   --> class: PC, n_carbon: 38, n_unsat: 3,
                                fa_comp: ((n_carbon: 18, n_unsat: 1), (n_carbon: 20, n_unsat: 2))
        Also, handles special fatty acid notations (modifiers) used for ceramides and
        plasmalogen lipids, examples:
            Cer(d36:2)      --> class: Cer, n_carbon: 36, n_unsat: 2, fa_mod: d
            Cer(d18:1/18:1) --> class: PC, n_carbon: 38, n_unsat: 3, fa_mod: d,
                                fa_comp: ((n_carbon: 18, n_unsat: 1), (n_carbon: 18, n_unsat: 1))
            PE(p40:4)       --> class: PE, n_carbon: 40, n_unsat: 4, fa_mod: p
            PE(p20:2/20:2)  --> class: PE, n_carbon: 40, n_unsat: 4, fa_mod: p,
                                fa_comp: ((n_carbon: 20, n_unsat: 2), (n_carbon: 20, n_unsat: 2))
        lipid name must conform to the general format:
            <lipid_class>([modifier]<n_carbon>:<n_unsat>[/<n_carbon>:<n_unsat>[/<n_carbon>:<n_unsat>]])
    parameters:
        name (str) -- lipid name to parse
    returns:
        (dict or None) -- parsed lipid information (always contains 'class', 'n_carbon', and 'n_unsat'
                    attributes) or None if it cannot be parsed as a lipid
"""
    parsed = {}

    # compile regex pattern
    l_pat = re.compile(
        r"^(?P<cls>[A-Za-z123]+)\((?P<mod>[pdoe]*)(?P<fc1>[0-9]+):(?P<fu1>[0-9]+)/*((?P<fc2>[0-9]+):(?P<fu2>[0-9]+))*/*((?P<fc3>[0-9]+):(?P<fu3>[0-9]+))*\)")

    # parse the name using regex
    l_res = l_pat.match(name)
    if l_res:
        # lipid class (required)
        if l_res.group('cls'):
            parsed["lipid_class"] = l_res.group('cls')
        else:
            # msg = "parse_lipid: failed to parse lipid class for: {}".format(name)
            # raise ValueError(msg)
            return None

        # value error due to failure to parse fatty acid composition
        # def raise_fatty_acid_value_error():
        #    msg = "parse_lipid: failed to parse fatty acid composition for: {}".format(name)
        #    raise ValueError(msg)

        # fc1 and fu1 are always required
        if not l_res.group('fc1') or not l_res.group('fu1'):
            # raise_fatty_acid_value_error()
            return None

        # check if a second fatty acid composition is supplied, e.g. (18:1/16:0)
        # if so, need to compute total fatty acid composition and add individual
        # fatty acids to a list
        if l_res.group('fc2'):
            if not l_res.group('fu2'):
                # raise_fatty_acid_value_error()
                return None
            # add info from the first two fatty acid compositions
            fc1, fu1 = int(l_res.group('fc1')), int(l_res.group('fu1'))
            fc2, fu2 = int(l_res.group('fc2')), int(l_res.group('fu2'))
            parsed["fa_comp"] = [
                {"n_carbon": fc1, "n_unsat": fu1},
                {"n_carbon": fc2, "n_unsat": fu2}
            ]
            # check for 3rd FA composition
            fc3, fu3 = 0, 0
            if l_res.group('fc3'):
                if not l_res.group('fu3'):
                    # raise_fatty_acid_value_error()
                    return None
                fc3, fu3 = int(l_res.group('fc3')), int(l_res.group('fu3'))
                parsed["fa_comp"].append({"n_carbon": fc3, "n_unsat": fu3})
            # compute total fatty acid composition
            parsed["n_carbon"] = fc1 + fc2 + fc3
            parsed["n_unsat"] = fu1 + fu2 + fc3
        else:
            # fc1 and fu1 are the total fatty acid composition
            parsed["n_carbon"] = int(l_res.group('fc1'))
            parsed["n_unsat"] = int(l_res.group('fu1'))

        # add fatty acid modifier if present
        if l_res.group('mod'):
            parsed["fa_mod"] = l_res.group('mod')
    else:
        # could not parse name as a lipid
        parsed = None

    return parsed


def fetch_lipid_class_log2fc(lipid_class, dataset, group_names, normed=False):
    """
fetch_lipid_class_log2fc
    description:
        Uses lipid identifications and already computed log2(fold-change) to gather all of the relevant data for a
        particular lipid class. Lipid identification must have already been performed and the group_names and normed
        args must be the same as the ones used to compute the original statistic.

        * when multiple MS adducts are present for a single lipid, the log2fas are averaged together *
    parameters:
        lipid_class (str) -- lipid class to analyze
        dataset (lipydomics.data.Dataset) -- lipidomics dataset
        group_names (list(str)) -- groups used to compute fold-change, only 2 groups allowed
        [normed (bool)] -- Use normalized data (True) or raw (False) [optional, default=False]
    returns:
        (tuple(list(), list(), list()) or tuple(None, None, None)) -- nc, nu, and log2fc for the specified lipid_class or None if no
                                                    corresponding data was found
"""
    # fetch the lipid annotations and log2fcs for lipids matching the lipid class
    log2fa_label = 'LOG2FC_{}_{}'.format('-'.join(group_names), 'normed' if normed else 'raw')
    lipids = {}
    for put_id, log2fa in zip(dataset.feat_ids, dataset.stats[log2fa_label]):
        if type(put_id) == list:  # actual identifications are lists not str
            parsed = parse_lipid(put_id[0].split('_')[0])
            if parsed is not None:  # putative id was parseable
                if parsed['lipid_class'] == lipid_class:  # lipid was the right lipid class
                    if log2fa not in [nan, inf]:  # log2fa cannot be nan or inf
                        nc, nu = parsed['n_carbon'], parsed['n_unsat']
                        if (nc, nu) not in lipids:
                            lipids[(nc, nu)] = [log2fa]
                        else:
                            lipids[(nc, nu)].append(log2fa)

    # average together multiple values, and split into separate lists
    nc, nu, l2 = [], [], []
    for comp in lipids:
        nc.append(comp[0])
        nu.append(comp[1])
        l2.append(mean(lipids[comp]) if len(lipids[comp]) > 1 else lipids[comp][0])

    return (None, None, None) if nc == [] else (nc, nu, l2)


def gen_tstamp():
    """
gen_tstamp
    description:
        generate a timestamp used for the build log and lipid database
    returns:
        (str) -- timestamp in YYMMDDhhmm format
"""
    now = datetime.now()
    s = '{:02d}{:02d}{:02d}{:02d}{:02d}'
    return s.format(now.year % 100, now.month, now.day, now.hour, now.minute)


def print_and_log(message, log_file, **kwargs):
    """
print_and_log
    description:
        prints a message to console and a specified log file
    parameters:
        message (str) -- message to print/write to log file
        log_file (file obj.) -- a file handle corresponding to the log file, should be opened in append mode
        **kwargs -- pass all other kwargs along to print(...)
"""
    print(message, **kwargs)
    print(message, file=log_file, **kwargs)


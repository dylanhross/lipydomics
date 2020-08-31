"""
    fill_measured_from_src.py
    Dylan H. Ross
    2018/10/04

        Fills the 'measured' table in lipids.db with values from the the cleaned reference datasets.
        After this script runs, the database will contain only relevant data on the measurements of these lipids, along
        with some metadata describining how they were measured
"""


import os
from sqlite3 import connect
from json import load as jload

from ..util import parse_lipid, print_and_log
from .build_params import include_ref_dsets


def add_src_dataset(cursor, src_tag, metadata, gid_start=0):
    """
add_src_dataset
    description:
        Adds values from a source dataset to the database, specified by a source tag
    parameters:
        cursor (sqlite3.cursor) -- cursor for running queries against the drugs.db database
        src_tag (str) -- source tag 
        metadata (dict(...)) -- CCS metadata: CCS type and method
        [gid_start (int)] -- starting number for m_id integer identifier [optional, default=0]
    returns:
        (int) -- the next available m_id value
"""
    ref_file = os.path.join(os.path.dirname(__file__), "reference_data/{}.json".format(src_tag))
    with open(ref_file, "r") as j:
        jdata = jload(j)
    # query string
    # m_id, name, adduct, mz, ccs, smi, src_tag
    qry = "INSERT INTO measured VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    # s_id starts at 0 and goes up from there
    m_id = gid_start
    for cmpd in jdata:

        # fix messed up adducts on the fly...
        adduct = cmpd["adduct"]
        fixed_adducts = {
            "[M+]+": "[M]+", 
            "M+NH4]+": "[M+NH4]+",
            "[M+H]+*": "[M+H]+",
            "[M+Na]+*": "[M+Na]+",
            "[M+H20-H]-": "[M+H2O-H]-"
        }
        if adduct in fixed_adducts:
            adduct = fixed_adducts[adduct]

        # use smi if included
        smi = None
        if "smi" in cmpd:
            smi = cmpd["smi"]
        
        # use rt if included
        rt = cmpd['rt'] if 'rt' in cmpd else None

        # add lipid class and FA composition information in by parsing the name
        parsed = parse_lipid(cmpd['name'])
        if not parsed:
            #print('name: {} from src_tag: {} not parsed as a lipid'.format(cmpd['name'], src_tag))
            # quitely skip unparsable lipids...
            pass
        else:
            l_cl, l_nc, l_nu = parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat']

            fa_mod = parsed['fa_mod'] if 'fa_mod' in parsed else None

            # CCS metadata
            ccs_type, ccs_method = metadata["type"], metadata["method"]

            qdata = (
                m_id, cmpd["name"], l_cl, l_nc, l_nu, fa_mod, adduct, cmpd["mz"], cmpd["ccs"],
                rt, smi, src_tag, ccs_type, ccs_method
            )
            cursor.execute(qry, qdata)
            m_id += 1

    return m_id


def main(tstamp):
    """ main build function """

    # connect to database
    db_path = os.path.join(os.path.dirname(__file__), 'lipids.db')
    con = connect(db_path)
    cur = con.cursor()

    # include reference datasets defined in build_params
    dsets = include_ref_dsets

    # CCS metadata by source
    metadata = {
        "zhou0817": {"type": "DT", "method": "single field, calibrated with Agilent tune mix (Agilent)"},
        "hine1217": {"type": "TW", "method": "calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-)"},
        "hine0217": {"type": "TW", "method": "calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-)"},
        "hine0119": {"type": "TW", "method": "calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-)"},
        'leap0219': {"type": "DT", "method": "stepped-field"},
        'blaz0818': {'type': 'DT', "method": 'single field'},
        'vasi0120_pos': {'type': 'TIMS', 'method': 'calibrated with 4 ions from ESI LC/MS tuning mix (Agilent)'},
        'vasi0120_neg': {'type': 'TIMS', 'method': 'calibrated with 4 ions from ESI LC/MS tuning mix (Agilent)'},
        'vasi0120_neg_corr': {'type': 'TIMS', 'method': 'calibrated with 4 ions from ESI LC/MS tuning mix (Agilent) *linear correction applied to CCS*'},
        'tsug0220_pos': {'type': 'TIMS', 'method': 'single field, calibrated'},
        'tsug0220_neg': {'type': 'TIMS', 'method': 'single field, calibrated'},
        'tsug0220_neg_corr': {'type': 'TIMS', 'method': 'single field, calibrated *linear correction applied to CCS'},
        'hine0520': {'type': 'TW', 'method': 'calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-)'}
    }

    # add each src dataset
    build_log = os.path.join(os.path.dirname(__file__), 'builds/build_log_{}.txt'.format(tstamp))
    with open(build_log, 'w') as bl:  # this is the first build script so OK to clear out the log file
        print_and_log("\nadding cleaned datasets into lipids.db", bl)

        gid_next = 0
        for dset in dsets:
            print_and_log("\tadding dataset: {} ...".format(dset), bl, end=" ")
            gid_next = add_src_dataset(cur, dset, metadata[dset], gid_start=gid_next)
            print_and_log("ok", bl)
        print_and_log("", bl)  # add a blank line

    # save changes to the database
    con.commit()
    con.close()

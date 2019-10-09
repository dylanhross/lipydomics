#!/usr/bin/python3
"""
    fill_theo_mz_from_gen.py
    Dylan H. Ross
    2019/10/07

        Fills the `theoretical_mz` table from lipids.db using enumeration over lipid classes, fatty acid composition, and
        MS adducts. 
"""


from generation import enumerate_all_lipids
from lipid_parser import parse_lipid


def add_enumerated_mz(cursor):
    """
add_src_dataset
    description:
        Adds enumerated m/z values into the database
    parameters:
        cursor (sqlite3.cursor) -- cursor for running queries against the drugs.db database
"""

    # query string
    # t_id, name, adduct, mz
    qry = 'INSERT INTO theoretical_mz VALUES (?,?,?,?,?,?,?,?)'
    # t_id starts at 0 and goes up from there
    t_id = 0
    for name, adduct, mz in enumerate_all_lipids():
        # parse the lipid name for class, nc, nu and FA mod
        parsed = parse_lipid(name)
        lc, nc, nu = parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat']
        fa_mod = parsed['fa_mod'] if 'fa_mod' in parsed else None
        qdata = (t_id, name, lc, nc, nu, fa_mod, adduct, mz)
        cursor.execute(qry, qdata)
        t_id += 1


if __name__ == '__main__':

    from sqlite3 import connect

    # connect to database
    con = connect("lipids.db")
    cur = con.cursor()

    print('adding theoretical m/z into lipds.db ...', end=' ')
    add_enumerated_mz(cur)
    print('ok')
    print()

    # save changes to the database
    con.commit()
    con.close()

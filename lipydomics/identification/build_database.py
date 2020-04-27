"""
    build_database.py
    Dylan H. Ross
    2020/02/29

        Replaces the old shell script used for building a new lipids database. Main method runs all of the individual
        build scripts.
"""


import os
import shutil
from sqlite3 import connect

from ..util import gen_tstamp
from .db_table_defs import measured, theo_mz, theo_ccs, theo_rt
from .fill_measured_from_src import main as fill_from_src
from .fill_theo_mz_from_gen import main as gen_theo_mz
from .train_lipid_ccs_pred import main as train_ccs_pred
from .characterize_lipid_ccs_pred import main as charac_ccs_pred
from .train_lipid_rt_pred import main as train_rt_pred
from .characterize_lipid_rt_pred import main as charac_rt_pred


def remove_old_files():
    """ removes old files that are going to be remade """

    # remove lipids.db
    this_dir = os.path.dirname(__file__)
    if os.path.isfile(os.path.join(this_dir, 'lipids.db')):
        os.remove(os.path.join(this_dir, 'lipids.db'))

    # remove all of the CCS and RT prediction performance images
    for d in ['ccs_pred_perf', 'rt_pred_perf']:
        img_dir = os.path.join(this_dir, d)
        for img in os.listdir(img_dir):
            os.remove(os.path.join(img_dir, img))


def initialize_db():
    """ initialize all of the tables in the database """

    db_file = os.path.join(os.path.dirname(__file__), 'lipids.db')
    con = connect(db_file)
    cur = con.cursor()
    for table in [measured, theo_mz, theo_ccs, theo_rt]:
        cur.execute(table)
    con.commit()
    con.close()


def make_database_copy(tstamp):
    """ makes a copy of the database and stores it in the builds directory """
    src = os.path.join(os.path.dirname(__file__), 'lipids.db')
    dst = os.path.join(os.path.dirname(__file__), 'builds/lipids_{}.db'.format(tstamp))
    shutil.copy(src, dst)


def main(tstamp):
    """ set up and run all of the individual build scripts """

    # get rid of old files
    remove_old_files()

    # initialize the database
    initialize_db()

    # run all of the build scripts
    fill_from_src(tstamp)
    gen_theo_mz(tstamp)
    train_ccs_pred(tstamp)
    charac_ccs_pred(tstamp)
    train_rt_pred(tstamp)
    charac_rt_pred(tstamp)

    # make a copy of the database and store it in the builds directory
    make_database_copy(tstamp)


# run the main function if this module is called directly
if __name__ == '__main__':
    main(gen_tstamp())


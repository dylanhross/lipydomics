#!/bin/bash
# create a new sqlite3 database using a SQL schema file then run the python initialization scripts
# to fill it with data from files in the reference_data directory and other generated data


# get rid of the database file (if it exists)
rm -f lipids.db

# initialize the database using the SQL schema file
cat lipids_schema.sql | sqlite3 lipids.db

# fill the database with measured values from the reference files
./fill_measured_from_src.py

# fill the theoretical mass table using enumeration
./fill_theo_mz_from_gen.py

# fill the theoretical CCS table using ML model trained on reference data
rm -f lipid_ccs_pred.pickle lipid_ccs_scale.pickle # get rid of pre-trained models (if present)
./train_lipid_ccs_pred.py


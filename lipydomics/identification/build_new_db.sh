#!/bin/bash
# create a new sqlite3 database using a SQL schema file then run the python initialization scripts
# to fill it with data from files in the cleaned_data directory and pull any missing information
# from other sources


# get rid of the database file (if it exists)
rm -f lipids.db

# initialize the database using the SQL schema file
cat lipids_schema.sql | sqlite3 lipids.db

# fill the database with initial values from the reference files
./fill_db_from_src.py


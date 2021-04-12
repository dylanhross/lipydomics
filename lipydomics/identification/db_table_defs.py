"""
    db_table_defs.py
    Dylan H. Ross
    2020/02/29

        CREATE TABLE statements for the lipids SQLite3 database
"""


measured = """
-- table in which to combine all of the measured datasets
CREATE TABLE measured (
    -- unique integer identifier
    m_id INTEGER UNIQUE NOT NULL,
    -- lipid name
    name TEXT NOT NULL,
    -- lipid_class and sum composition
    lipid_class TEXT NOT NULL,
    lipid_nc INTEGER NOT NULL,
    lipid_nu INTEGER NOT NULL,
    -- fatty acid modifier
    fa_mod TEXT,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- m/z and CCS
    mz REAL NOT NULL,
    ccs REAL NOT NULL,
    -- HILIC rt
    rt REAL,
    -- neutral smiles structure
    smi TEXT,
    -- tag referencing which dataset the value is from
    src_tag TEXT NOT NULL,
    -- CCS type (DT, TW, ...)
    ccs_type TEXT NOT NULL,
    -- describe method used for CCS measurement (e.g. stepped-field, calibrated with polyalanine)
    ccs_method TEXT
);
"""

theo_mz = """
-- table to store predicted masses from enumeration on composition and MS adducts
CREATE TABLE predicted_mz (
    -- unique integer identifier
    t_id INTEGER UNIQUE NOT NULL,
    -- lipid name
    name TEXT NOT NULL,
    -- lipid_class and sum composition
    lipid_class TEXT NOT NULL,
    lipid_nc INTEGER NOT NULL,
    lipid_nu INTEGER NOT NULL,
    -- fatty acid modifier
    fa_mod TEXT,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- predicted m/z
    mz REAL NOT NULL
);
"""

theo_ccs = """
-- table to store predicted CCS values, associated with predicted m/z values
CREATE TABLE predicted_ccs (
    -- unique integer identifier
    t_id INTEGER UNIQUE NOT NULL,
    -- predicted CCS
    ccs REAL NOT NULL
);
"""

theo_rt = """
-- table to store predicted retention time values, associated with predicted m/z values
CREATE TABLE predicted_rt (
    -- unique integer identifier
    t_id INTEGER UNIQUE NOT NULL,
    -- predicted retention time
    rt REAL NOT NULL
);
"""


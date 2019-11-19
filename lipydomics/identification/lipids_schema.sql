

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


-- table to store theoretical masses from enumeration on composition and MS adducts
CREATE TABLE theoretical_mz (
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
    -- theoretical m/z
    mz REAL NOT NULL
);


-- table to store theoretical CCS values, associated with theoretical m/z values
CREATE TABLE theoretical_ccs (
    -- unique integer identifier
    t_id INTEGER UNIQUE NOT NULL,
    -- theoretical m/z
    ccs REAL NOT NULL
);


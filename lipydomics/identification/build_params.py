"""
    build_params.py
    Dylan H. Ross
    2020/03/11

    description:
        defines parameters to use for building the database:
            include_ref_dsets - reference datasets to include in the database
            ccs_pred_ref_dsets - reference datasets to use for CCS prediction
"""


# reference datasets to include in the 'measured' table
include_ref_dsets = [
    "zhou0817",
    "hine1217",
    "hine0217",
    "hine0119",
    'leap0219',
    'blaz0818',
    'vasi0120_pos',
    'vasi0120_neg',
    'vasi0120_neg_corr',
    'tsug0220_pos',
    'tsug0220_neg',
    'tsug0220_neg_corr',
    'hine0520'
]


# reference datasets to use for CCS prediction
ccs_pred_ref_dsets = [
    "zhou0817",
    "hine1217",
    "hine0217",
    "hine0119",
    'leap0219',
    'blaz0818',
    'vasi0120_pos',
    'vasi0120_neg_corr'
    'tsug0220_pos',
    'tsug0220_neg_corr',
    'hine0520'
]


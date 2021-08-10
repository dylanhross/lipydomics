"""
    build_params.py
    Dylan H. Ross
    2020/03/11

    description:
        defines parameters to use for building the database:
            include_ref_dsets - reference datasets to include in the database
            ccs_pred_ref_dsets - reference datasets to use for CCS prediction


    NOTE:
        Vasilopoulou et al. data excluded from database build and CCS predictive model
        based on issues raised in:

        Nat. Comm. (2021) 12:4771
        Quality control requirements for the correctannotation of lipidomics data
        Harald C. Köfeler, Thomas O. Eichmann, Robert Ahrends, John A. Bowden,
        Niklas Danne-Rasche, Edward A. Dennis, Maria Fedorova, William J. Griffiths,
        Xianlin Han, Jürgen Hartler, Michal Holčapek, Robert Jirásko, Jeremy P.
        Koelmel, Christer S. Ejsing, Gerhard Liebisch, Zhixu Ni,
        Valerie B. O’Donnell, Oswald Quehenberger, Dominik Schwudke
        Andrej Shevchenko, Michael J. O. Wakelam, Markus R. Wenk, Denise Wolrab
        & Kim Ekroos

        ARISING  FROM M. Mann et al. Nature  Communications https://doi.org/10.1038/s41467-019-14044-x(2020).
"""


# reference datasets to include in the 'measured' table
include_ref_dsets = [
    "zhou0817",
    "hine1217",
    "hine0217",
    "hine0119",
    'leap0219',
    'blaz0818',
    # See note in description...
    #'vasi0120_pos',
    #'vasi0120_neg',
    #'vasi0120_neg_corr',
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
    # See note in description...
    #'vasi0120_pos',
    #'vasi0120_neg_corr'
    'tsug0220_pos',
    'tsug0220_neg_corr',
    'hine0520'
]


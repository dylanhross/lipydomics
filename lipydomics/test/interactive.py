"""
    lipydomics/test/interactive.py
    Dylan H. Ross
    2020/01/17

    description:
        Define tests associated with the lipydomics/interactive.py module
"""


from lipydomics.interactive import abbreviate_sheet


def abbrev_xl_sheet_names():
    """
abbrev_xl_sheet_names
    description:
        Tests the function that abbreviates sheet names for exporting a Dataset into an Excel Workbook to ensure that
        they are <= 31 characters in length

        Test fails if there are any errors, or if any of the abbreviated names are > 31 characters in length
    returns:
        (bool) -- test pass (True) or fail (False)
"""
    names_to_test = [
        '',  # empty
        'PCA3_A-B-C-D-E_projections_raw',  # already short enough
        'PCA3_A-B-C-D-E-F-G-H_projections_raw',  # needs to be abbreviated, group names already short
        'PCA3_A-B-C-D-E-F-G-H-I-J-K-L-M-N_projections_raw',  # same, but too long so ends up just truncated
        'PCA3_Alabama-Bigbird-Cauliflower_loadings_normed',  # needs to be abbreviated, group names as well
        'PLS-DA_Bigbird-Cauliflower_loadings_normed',  # same but different stat
        'ANOVA_Alabama-Bigbird-Cauliflower-Dingle_raw',  # same but different stat
        '2-group-corr_Bigbird-Cauliflower_normed',  # same but different stat
        'ThisShouldBeTruncatedToOnly31Characters'  # unrecognized nonsense, should just be truncated
    ]
    expected_abbreviations = [
        '',  # empty, nothing to do
        'PCA3_A-B-C-D-E_projections_raw',  # already short enough, nothing to do
        'PC3PR_A-B-C-D-E-F-G-H',  # abbreviated, group names already short
        'PCA3_A-B-C-D-E-F-G-H-I-J-K-L-M-',  # same, but too long so ends up just truncated
        'PC3LN_Ala-Big-Cau',  # abbreviated, group names too
        'PLSLN_Big-Cau',  # same but different stat
        'ANVR_Ala-Big-Cau-Din',  # same but different stat
        '2GCN_Big-Cau',  # same but different stat
        'ThisShouldBeTruncatedToOnly31Ch'  # truncate to 31 characters
    ]
    actual_abbreviations = [abbreviate_sheet(name) for name in names_to_test]
    for act, exp in zip(actual_abbreviations, expected_abbreviations):
        if act != exp:
            print('abbrev_xl_sheet_names\n\texpected: "{}"\n\tactual: "{}"\n'.format(exp, act))
            return False
    return True

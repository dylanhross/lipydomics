"""
    lipydomics/test/__main__.py
    Dylan H. Ross
    2019/02/02

    description:
        Calls the run_all_tests method to run all of the tests defined in the various testing modules
"""


import traceback
import os




from lipydomics.test.data import dataset_init_mock1, dataset_normalize_mock1, dataset_getgroup_mock1
from lipydomics.test.stats import addanovap_mock1, addpca3_mock1, addplsda_mock1, addplsda_3groups_mock1
from lipydomics.test.plotting import (
    barplot_feature_bygroup_mock1, scatter_pca3_projections_bygroup_mock1, scatter_plsda_projections_bygroup_mock1,
    splot_plsda_pcorr_bygroup_mock1
)


def run_all_tests():
    """
run_all_tests
    description:
        runs all tests sequentially, if there are any failures the test function docstring is printed. A traceback is 
        printed as well if an exception is thrown.
"""
    # references to al of the test functions to be run, and order to run them in
    all_tests = [
        dataset_init_mock1,
        dataset_normalize_mock1,
        dataset_getgroup_mock1,
        addanovap_mock1,
        addpca3_mock1,
        addplsda_mock1,
        addplsda_3groups_mock1,
        barplot_feature_bygroup_mock1,
        scatter_pca3_projections_bygroup_mock1,
        scatter_plsda_projections_bygroup_mock1,
        splot_plsda_pcorr_bygroup_mock1
    ]
    # run the tests
    failed = False
    print("running all tests ... ", end="")
    for test in all_tests:
        try:
            passed = test()
        except Exception as e:
            print('\n', test.__doc__, 'TEST FAILED WITH EXEPTION!\n', e)
            print(traceback.format_exc())
            failed = True
            break
        if not passed:
            print('\n', test.__doc__, 'TEST FAILED!\n')
            failed = True
            break
    if not failed:
        print("passed")



if __name__ == '__main__':
    run_all_tests()

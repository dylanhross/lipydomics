"""
    lipydomics/test/__init__.py
    Dylan H. Ross
    2019/02/02

    description:
        TODO
"""

import traceback


def run_tests(tests):
    """
run_tests
    description:
        runs all specified test functions sequentially, if there are any failures the test function docstring is
        printed. A traceback is printed as well if an exception is thrown.
    paramters:
        tests (list(function)) -- a list of all individual test functions to run
"""
    # run the tests
    failed = False
    print("running all tests ... ", end="")
    for test in tests:
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


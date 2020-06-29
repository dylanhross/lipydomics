"""
    lipydomics/test/__main__.py
    Dylan H. Ross
    2019/02/02

    description:
        Calls the run_all_tests method to run all of the tests defined in the various testing modules
"""


from lipydomics.test import run_tests
from lipydomics.test.data import all_tests as data_all_tests
from lipydomics.test.stats import all_tests as stats_all_tests
from lipydomics.test.plotting import all_tests as plotting_all_tests
from lipydomics.test.identification import all_tests as identification_all_tests
from lipydomics.test.rt_calibration import all_tests as rt_calibration_all_tests
from lipydomics.test.util import all_tests as util_all_tests


# references to al of the test functions to be run, and order to run them in
all_tests = data_all_tests + stats_all_tests + plotting_all_tests + identification_all_tests + \
            rt_calibration_all_tests + util_all_tests
run_tests(all_tests)

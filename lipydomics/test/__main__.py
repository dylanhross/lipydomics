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


# run through each module's tests
print('(lipydomics.test.data) ', end='', flush=True)
run_tests(data_all_tests)

print('(lipydomics.test.stats) ', end='', flush=True)
run_tests(stats_all_tests)

print('(lipydomics.test.plotting) ', end='', flush=True)
run_tests(plotting_all_tests)

print('(lipydomics.test.identification) ', end='', flush=True)
run_tests(identification_all_tests)

print('(lipydomics.test.rt_calibration) ', end='', flush=True)
run_tests(rt_calibration_all_tests)

print('(lipydomics.test.util) ', end='', flush=True)
run_tests(util_all_tests)

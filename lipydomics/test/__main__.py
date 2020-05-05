"""
    lipydomics/test/__main__.py
    Dylan H. Ross
    2019/02/02

    description:
        Calls the run_all_tests method to run all of the tests defined in the various testing modules
"""


import traceback


from lipydomics.test.data import (
    dataset_init_mock1, dataset_normalize_mock1, dataset_getgroup_mock1, dataset_assign_groups_using_replicates_real1,
    dataset_save_load_bin_mock1, dataset_export_feature_data_real1, dataset_export_xlsx_real1,
    dataset_export_analyzed_xlsx_real1
)
from lipydomics.test.stats import (
    addanovap_mock1, addanovap_real1, addpca3_mock1, addpca3_real1, addplsda_mock1, addplsda_3groups_mock1, 
    addplsda_real1, add2groupcorr_mock1, add2groupcorr_3groups_mock1, add2groupcorr_real1, addplsra_real1,
    addlog2fc_real1
)
from lipydomics.test.plotting import (
    barplot_feature_bygroup_mock1, batch_barplot_feature_bygroup_real1, scatter_pca3_projections_bygroup_mock1,
    scatter_plsda_projections_bygroup_mock1, splot_plsda_pcorr_bygroup_mock1, scatter_plsra_projections_bygroup_real1,
    heatmap_lipid_class_log2fc_real1
)
from lipydomics.test.identification import (
    add_feature_ids_any_real1, add_feature_ids_any_real1_tstamp, add_feature_ids_real1_bad_tstamp, predict_ccs_noerrs,
    predict_ccs_notencodable, predict_rt_noerrs, predict_rt_notencodable
)
from lipydomics.test.rt_calibration import (
    get_ref_rt_lipids1, rtcal_init_mismatch_len, rtcal_calibrate_rtc1_c12, rtcal_calibrate_rtc1_c13,
    rtcal_calibrate_rtc2_c12, rtcal_calibrate_rtc2_c13
)
from lipydomics.test.util import abbrev_xl_sheet_names, fetch_lipid_class_log2fa_real1


def run_all_tests():
    """
run_all_tests
    description:
        runs all tests sequentially, if there are any failures the test function docstring is printed. A traceback is 
        printed as well if an exception is thrown.
"""
    # references to al of the test functions to be run, and order to run them in
    all_tests = [
        # test/data
        dataset_init_mock1,
        dataset_normalize_mock1,
        dataset_getgroup_mock1,
        dataset_assign_groups_using_replicates_real1,
        dataset_save_load_bin_mock1,
        dataset_export_feature_data_real1,
        dataset_export_xlsx_real1,
        dataset_export_analyzed_xlsx_real1,
        # test/stats
        addanovap_mock1,
        addanovap_real1,
        addpca3_mock1,
        addpca3_real1,
        addplsda_mock1,
        addplsda_3groups_mock1,
        addplsda_real1,
        add2groupcorr_mock1,
        add2groupcorr_3groups_mock1,
        add2groupcorr_real1,
        addplsra_real1,
        addlog2fc_real1,
        # test/plotting
        barplot_feature_bygroup_mock1,
        batch_barplot_feature_bygroup_real1,
        scatter_pca3_projections_bygroup_mock1,
        scatter_plsda_projections_bygroup_mock1,
        splot_plsda_pcorr_bygroup_mock1,
        scatter_plsra_projections_bygroup_real1,
        heatmap_lipid_class_log2fc_real1,
        # test/identification
        add_feature_ids_any_real1,
        add_feature_ids_any_real1_tstamp,
        add_feature_ids_real1_bad_tstamp,
        predict_ccs_noerrs,
        predict_ccs_notencodable,
        predict_rt_noerrs,
        predict_rt_notencodable,
        # test/rt_calibration
        get_ref_rt_lipids1,
        rtcal_init_mismatch_len,
        rtcal_calibrate_rtc1_c12,
        rtcal_calibrate_rtc1_c13,
        rtcal_calibrate_rtc2_c12,
        rtcal_calibrate_rtc2_c13,
        # test/util
        abbrev_xl_sheet_names,
        fetch_lipid_class_log2fa_real1
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

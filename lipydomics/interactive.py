"""
    lipydomics/interactive.py
    Jang Ho Cho and Dylan Ross

    description:

"""


from pandas import DataFrame, concat, read_csv
import csv
import numpy as np
import os

from lipydomics.data import Dataset
from lipydomics.stats import add_anova_p, add_pca3, add_plsda, add_2group_corr, add_plsra, add_log2fc
from lipydomics.plotting import (
    barplot_feature_bygroup, batch_barplot_feature_bygroup, scatter_pca3_projections_bygroup,
    scatter_plsda_projections_bygroup, splot_plsda_pcorr_bygroup, scatter_plsra_projections_bygroup,
    heatmap_lipid_class_log2fc
)
from lipydomics.identification import add_feature_ids, remove_potential_nonlipids
from lipydomics.identification.rt_calibration import get_ref_rt, RTCalibration
from lipydomics.util import filter_d, parse_lipid


def load_dset():
    """
load_dset
    description:
        Prompts the user with options to load a lipydomics Dataset instance, either a new one from a .csv file or an
        existing one from .pickle file. Returns the loaded instance or None on any sort of failure. There is also a
        hidden exit option that will completely stop execution.
        (if None is returned, this function is called again until a Dataset instance is returned)
        (if 'exit' is returned, main() will return None rather than actually exiting)
    returns:
        (lipydomics.data.Dataset or None or 'exit') -- lipidomics dataset instance
"""
    dset = None
    print('\nWhat would you like to do?')
    print('\t1. Make a new Dataset')
    print('\t2. Load a previous Dataset')
    option = input('> ')

    if option == '1':
        print('Please enter the path to the csv file you want to work with.')
        csv_fname = input('> ')
        # validate the csv file exists
        if not os.path.isfile(csv_fname):
            print('! ERROR: Make sure the path specified is correct and the file exists.')
            return None
        # prompt for positive or negative ESI mode
        print('What ESI mode was used for this data? (pos/neg)')
        esi = input('> ')
        if esi not in ['pos', 'neg']:
            print('! ERROR: ESI mode "{}" not recognized'.format(esi))
            return None
        # load the Dataset from .csv file
        dset = Dataset(csv_fname, esi_mode=esi)
        print('! INFO: Loaded a new Dataset from .csv file: "{}"'.format(csv_fname))
        # try to automatically assign headers
        print('Would you like to automatically assign groups from headers? (y/N)')
        ans = input('> ')
        if ans == 'y':
            try:
                with open(csv_fname, newline='') as f:
                    reader = csv.reader(f)
                    header = next(reader)
                header = header[3:]
                group_map = {}
                for i in range(len(header)):
                    if header[i] not in group_map:
                        group_map[header[i]] = [i]
                    else:
                        group_map[header[i]] = group_map[header[i]] + [i]
                dset.assign_groups(group_map)
                print('! INFO: Automatically assigned groups from headers')
            except Exception as e:
                print('! ERROR:', e)
                print('! ERROR: Unable to automatically assign groups from headers')
                # reload the Dataset just in case
                dset = Dataset(csv_fname)

    elif option == '2':
        print('Please enter the path to the pickle file you want to load.')
        pickle_fname = input('> ')
        if not os.path.isfile(pickle_fname):
            print('! ERROR: Make sure the path specified is correct and the file exists.')
            return None
        dset = Dataset.load_bin(pickle_fname)
        print('! INFO: Loaded existing Dataset from .pickle file: "{}"'.format(pickle_fname))

    # exit option not listed
    elif option == 'exit':
        return 'exit'

    return dset


def manage_groups(dset):
    """
manage_groups
    description:
        Prompts the user with options to manage group assignments on the Dataset instance:
            - Assign indices to a group
            - View assigned group indices
            - Get data by group
        Returns a boolean indicating whether the user is finished with assigning groups.
        (this function gets called again if False is returned)
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- finished managing groups
"""
    print('Managing groups... What would you like to do?')
    print("\t1. Assign group")
    print("\t2. View assigned groups")
    print("\t3. Get data by group(s)")
    print('\t"back" to go back')
    option = input('> ')

    if option == "1":
        print("Please provide a name for a group and its indices in order of name > starting index > ending index."
              "\n\t* group name should not contain spaces\n\t* indices start at 0\n\t* example: 'A 1 3'")
        group = input('> ')
        group = group.split()
        name = group[0]
        indices = [_ for _ in range(int(group[1]), int(group[2]) + 1)]
        try:
            dset.assign_groups({name: indices})
            print('! INFO: Assigned indices: {} to group: "{}"'.format(dset.group_indices[name], name))
        except ValueError as ve:
            print('! ERROR:', ve)
            print("! ERROR: Failed to assign group, please check your formatting and try again")

    elif option == "2":
        for group in dset.group_indices:
            print('\t"{}": {}'.format(group, dset.group_indices[group]))
        return False

    elif option == "3":
        print("Which group would you like to view?")
        name = input('> ')
        print(dset.get_data_bygroup(name))
        return False

    elif option == 'back':
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


def filter_data(dset):
    """
filter_data
    description:
        Prompts the user with options for filtering the data:
            - The user is to provide ranges for M/Z, RT and CCS values and it let's the user download
              a csv file containing all data matching that range.
            - Can also take multiple ranges given a CSV file of ranges.
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- finished filtering data
"""
    print('Filtering data... What would you like to do?')
    print("\t1. Single query")
    print("\t2. Batch query")
    print("\t3. S-Plot filtering")
    print('\t"back" to go back')
    option = input('> ')
    if option == "back":
        return True
    label_dat = dset.labels
    label_df = DataFrame(label_dat)
    if option == "1":
        print("Please Provide m/z and tolerance (Ex. '150 1'  <--- This would be 150 plus or minus 1)")
        mz = input('> ')
        print("Please Provide Retention Time and tolerance (Ex. '1 1' <--- This would be 1 plus or minus 1)")
        rt = input('> ')
        print("Please Provide CCS and tolerance (Ex. '150 3'  <--- This would be 150 plus or minus 3%)")
        ccs = input('> ')
        print("Which group would you like to choose? ('All' to select the whole data)")
        group = input('> ')
        try:
            if group == "All":
                cur_data = DataFrame(dset.intensities)
            else:
                cur_data = dset.get_data_bygroup(group)
            int_df = DataFrame(cur_data)
            cur_df = concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
            mzs = [float(_) for _ in mz.split()]
            rts = [float(_) for _ in rt.split()]
            ccss = [float(_) for _ in ccs.split()]
            # convert CCS tolerance from percentage to an absolute value
            ccss[1] = ccss[1] / 100. * ccss[0]
            filtered = filter_d(mzs, rts, ccss, cur_df)
        except ValueError as ve:
            print('! ERROR:', ve)
            print("! ERROR: Failed to filter data, please check your groups and try again")
            return False

    elif option == "2":
        print("Please provide the path of the file with batch-query information")
        path = input('> ')
        try:
            query = read_csv(path)
        except Exception as e:
            print('! ERROR:', e)
            print("! ERROR: Failed to load the file. Please make sure the file exists at the right path.")
            return False

        print("Which group would you like to choose? ('All' to select the whole data)")
        group = input('> ')
        try:
            if group == "All":
                cur_data = DataFrame(dset.intensities)
                cur_df = concat([label_df, cur_data], axis=1, ignore_index=True, sort=False)
            else:
                cur_data = dset.get_data_bygroup(group)
                int_df = DataFrame(cur_data)
                cur_df = concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
        except ValueError as ve:
            print('! ERROR:', ve)
            print("! ERROR: Failed to filter data, please check your groups and try again")
            return False
        for index, row in query.iterrows():
            if index == 0:
                filtered = filter_d([row["m/z"], row["m/z_tol"]], [row["rt"], row["rt_tol"]],
                                    [row["ccs"], row["ccs_tol"]], cur_df)
            else:
                filtered = pd.concat([filtered,
                                      filter_d([row["m/z"], row["m/z_tol"]],
                                               [row["rt"], row["rt_tol"]], [row["ccs"], row["ccs_tol"]],
                                               cur_df)])
    elif option == "3":
        print("Which groups would you like to choose?")
        groups = input("> ")
        group_names = groups.split()
        print("Would you like to filter on normalized data? (y/N)")
        nrm = input("> ")
        if nrm == "y":
            nrm = "norm"
        else:
            nrm = "raw"
        try:
            x = dset.stats['PLS-DA_{}_loadings_{}'.format('-'.join(group_names), nrm)].T[0]
            y = dset.stats['2-group-corr_{}_{}'.format('-'.join(group_names), nrm)]
            x = abs(x)
            y = abs(y)
        except KeyError as ke:
            print('! ERROR:', ke)
            print("! ERROR: Required Statistics not yet computed")
            return False
        print("Please provide desired range on PLS-DA loadings separated by space. "
              "\n\tExample: '-10 10' will get values between -10 and 10")
        x_range = input("> ")
        x_range = [float(i) for i in x_range.split()]
        print("Please provide desired range on group correlation values. "
              "\n\tExample: '-10 10' will get values between -10 and 10")
        y_range = input("> ")
        y_range = [float(i) for i in y_range.split()]
        indices = []
        for i in range(0, len(x)):
            if x[i] <= x_range[1] and x[i] >= x_range[0] and y[i] <= y_range[1] and y[i] >= y_range[0]:
                indices.append(i)
        if not indices:
            print("! ERROR: Could not find any matching point")
        else:
            filtered = []
            group_data = list(dset.get_data_bygroup(group_names))
            first = group_data[0]
            second = group_data[1]
            for ind in indices:
                t = list(dset.labels[ind]) + list(first[ind]) + list(second[ind])
                filtered.append(t)
            filtered = DataFrame(filtered)

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False
    print("! INFO: Successfully filtered data. Would you like to download the result as csv? (y/N)")
    option = input('> ')
    if option == "y":
        print("Please specify the path you want to save the csv file")
        path = input('> ')
        try:
            filtered.to_csv(path, index=None, header=False)
            print("! INFO: Successfully downloaded the filtered data")
        except Exception as e:
            print('! ERROR:', e)
            print("! ERROR: Failed to download, please check your path and try again")
            return False

    return True


def manage_statistics(dset):
    """
manage_statistics
    description:
        Prompts the user with options to manage statistical analyses on the Dataset instance:
            -
        Returns a boolean indicating whether the user is finished with assigning groups.
        (this function gets called again if False is returned)
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- finished managing statistics
"""
    print('Managing statistics... What would you like to do?')
    print("\t1. Compute Statistics")
    print("\t2. View Statistics")
    print("\t3. Export .csv File of Computed Statistics")
    print('\t"back" to go back')
    option = input('> ')

    if option == '1':
        print('Computing statistics... What would you like to do?')
        print("\t1. Anova-P")
        print("\t2. PCA3")
        print("\t3. PLS-DA")
        print("\t4. Two Group Correlation")
        print("\t5. PLS-RA (using external continuous variable)")
        print("\t6. Two Group Log2(fold-change)")
        print('\t"back" to go back')
        opt2 = input('> ')
        # map options to stats functions
        stats_f_map = {
            '1': add_anova_p, '2': add_pca3, '3': add_plsda, '4': add_2group_corr, '5': add_plsra, '6': add_log2fc
        }
        if opt2 in stats_f_map:
            print('Would you like to use normalized data? (y/N)')
            norm = input('> ') == 'y'
            if norm:
                # make sure there is actually normalized data in the Dataset
                if dset.normed_intensities is None:
                    print('! ERROR: No normalization has been performed yet')
                    return False
            print('Please enter group names to use in this analysis, separated by spaces')
            groups = input('> ').split()
            try:
                if opt2 == '5':
                    # get the y variable
                    print('Please enter the path to a text file with the target variable values for regression')
                    print('\tExample: external.txt')
                    target_path = input('> ')
                    with open(target_path, 'r') as f:
                        y = np.array([float(_.strip()) for _ in f.readlines()])
                    stats_f_map[opt2](dset, groups, y, normed=norm)
                else:
                    stats_f_map[opt2](dset, groups, normed=norm)
                print('! INFO: Applied new statistical analysis using groups: {}'.format(groups))
            except ValueError as ve:
                print('! ERROR:', ve)
                print('! ERROR: Unable to perform statistical analysis, check group names and try again')
            except FileNotFoundError as fnfe:
                print('! ERROR:', fnfe)
                print('! ERROR: Unable to find file with target variable values: {}'.format(target_path))
        elif opt2 == 'back':
            pass
        else:
            print('! ERROR: unrecognized option: "{}"'.format(opt2))
        return False

    elif option == '2':
        if dset.stats:
            print("Computed statistics are:")
            for s in dset.stats:
                print("\t" + s + "\n")
        else:
            print("No statistic have been computed yet")
        return False

    elif option == '3':
        print('Please enter the path to save the exported statics file (.csv) under.')
        export_path = input('> ')
        try:
            pd.DataFrame.from_dict(dset.stats).to_csv(export_path)
            print('! INFO: Exported statistic to file: "{}"'.format(export_path))
        except Exception as e:
            print('! ERROR:', e)
            print('! ERROR: Failed to export statistics to file: "{}"'.format(export_path))
        return False

    elif option == 'back':
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


def make_plots(dset):
    """
make_plots
    description:

    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- finished making plots
"""
    print('Making Plots... What would you like to do?')
    print("\t1. Bar plot feature by group")
    print("\t2. Batch bar plot features by group")
    print("\t3. Scatter PCA3 Projections by group")
    print("\t4. Scatter PLS-DA Projections by group")
    print("\t5. S-Plot PLS-DA and Pearson correlation by group")
    print("\t6. Scatter PLS-RA Projections by group")
    print("\t7. Heatmap of Log2(fold-change) by lipid class")
    print('\t"back" to go back')
    option = input('> ')

    # map the options to plotting functions
    plots_f_map = {
        '1': barplot_feature_bygroup, '2': batch_barplot_feature_bygroup, '3': scatter_pca3_projections_bygroup,
        '4': scatter_plsda_projections_bygroup, '5': splot_plsda_pcorr_bygroup, '6': scatter_plsra_projections_bygroup,
        '7': heatmap_lipid_class_log2fc
    }
    if option in plots_f_map:
        print("Where would you like to save the plot(s)? (default = current directory)")
        plot_dir = input('> ')
        print("Which groups would you like to plot (separated by spaces)?")
        groups = input('> ').split()
        print('Would you like to use normalized data? (y/N)')
        norm = input('> ') == 'y'
        if norm:
            # make sure there is actually normalized data in the Dataset
            if dset.normed_intensities is None:
                print('! ERROR: No normalization has been performed yet')
                return False
        if option in ['1', '2']:
            try:
                if option == '1':
                    # get the m/z, rt, and CCS of the feature
                    print("Please enter the m/z, retention time, and CCS of the feature\n\t* separated by spaces\n\t*"
                          " example: '345.6729 1.23 123.45'")
                    feature = input('> ').split()
                    feat = [float(f) for f in feature]
                else:
                    # get the name of the input file
                    print('Please enter the path to the .csv file containing features to plot...',
                          "\n\t* example: 'plot_these_features.csv'")
                    feat = input('> ')  # store the .csv file name in feat
                    # check that the file exists
                    if not os.path.isfile(feat):
                        print('! ERROR: cannot find input file "{}", check path and try again'.format(feat))
                        return False
                # get the m/z, rt, and CCS tolerance
                print("Please enter the search tolerances for m/z, retention time, and CCS\n\t* separated by spaces"
                      "\n\t* example: '0.01 0.5 3.0'\n\t* CCS tolerance is a percentage, not an absolute value")
                tolerance = input('> ').split()
                tol = [float(t) for t in tolerance]
                result = plots_f_map[option](dset, groups, plot_dir, feat, normed=norm, tolerance=tol)
                if result or option == '2':
                    print('! INFO: Generated plot(s) for groups: {}'.format(groups))
                else:
                    print('! ERROR: Could not match feature')
            except ValueError as ve:
                print('! ERROR:', ve)
                print('! ERROR: Unable to generate plot using groups: {}'.format(groups))
        elif option == '7':
            try:
                # prompt for the lipid class
                print("Please enter the lipid class you would like to generate a heatmap with")
                lipid_class = input('> ')
                found_lipid = plots_f_map[option](lipid_class, dset, groups, plot_dir, normed=norm)
                if not found_lipid:
                    m = '! ERROR: unable to find lipids with matching identifications (lipid class: {})'
                    print(m.format(lipid_class))
                else:
                    print('! INFO: Generated heatmap for lipid class: {}'.format(lipid_class))
            except ValueError as ve:
                print('! ERROR:', ve)
                print('! ERROR: Unable to generate plot')
        else:
            try:
                plots_f_map[option](dset, groups, plot_dir, normed=norm)
                print('! INFO: Generated plot for groups: {}'.format(groups))
            except KeyError as ke:
                print('! ERROR:', ke)
                print('! ERROR: Required statistic(s) have not been computed for groups: {}'.format(groups))
            except ValueError as ve:
                print('! ERROR:', ve)
                print('! ERROR: Unable to generate plot using groups: {}'.format(groups))
        return False

    elif option == 'back':
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


def identify_lipids(dset):
    """
identify_lipids
    description:
        Prompts the user with options to perform lipid identification at a variety of levels of confidence
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
"""
    print("Identifying Lipids... Please enter the tolerances for m/z, retention time and CCS matching\n\t* separated by"
          " spaces\n\t* example: '0.01 0.5 3.0'\n\t* CCS tolerance is a percentage, not an absolute value")
    tolerance = input('> ').split()
    tol = [float(t) for t in tolerance]
    print('What type of tolerance should be used for m/z (Da/ppm) ?')
    mz_tol_type = input('> ')
    print("Please specify an identification level")
    print("\t'pred_mz' - match on predicted m/z")
    print("\t'pred_mz_rt' - match on predicted m/z and retention time")
    print("\t'pred_mz_ccs' - match on predicted m/z and CCS")
    print("\t'pred_mz_rt_ccs' - match on predicted m/z, retention time, and CCS")
    print("\t'meas_mz_ccs' - match on measured m/z and CCS")
    print("\t'meas_mz_rt_ccs' - match on measured m/z, retention time, and CCS")
    print("\t'any' - try all criteria (highest confidence first)")
    print("\t'back' to go back")
    option = input('> ')

    if option in ['pred_mz', 'pred_mz_rt', 'pred_mz_ccs', 'pred_mz_rt_ccs', 'meas_mz_ccs', 'meas_mz_rt_ccs', 'any']:
        # make the identifications
        n_identified = None
        try:
            add_feature_ids(dset, tol, level=option, mz_tol_type=mz_tol_type)
            n_identified = len([_ for _ in dset.feat_ids if type(_) == list])
            print("! INFO: Lipid identifications added successfully ({} lipids identified)".format(n_identified))
        except ValueError as ve:
            print('! ERROR:', ve)
            print("! ERROR: Unable to make lipid identifications")

        # optionally filter out dubious identifications after identifying
        if n_identified is not None and option in ['pred_mz', 'pred_mz_rt', 'any']:
            print('Would you like to remove identifications that do not look like lipids (based on CCS vs. m/z trend)?')
            option2 = input('> ')
            if option2 in ['y', 'Y', 'yes', 'YES', 'Yes']:
                n = remove_potential_nonlipids(dset)
                print('! INFO: removed {} identifications with CCS beyond +/-10% from global lipid trend'.format(n))

    elif option == 'back':
        # return None to finish
        return

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        # return None to finish
        return


def normalize_data(dset, df):
    """
normalize_data
    description:
        Prompts the user with options for normalizing the data:
            - The user is to choose between internal and external normalization and
              it will normalize the data within lipidomics dataset
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
        df (Pandas DataFrame) -- DataFrame version of lipidomics dataset
    returns:
        (bool) -- finished normalizing data
"""
    print('Normalizing data... What would you like to do?')
    print("\t1. Internal")
    print("\t2. External")
    print('\t"back" to go back')
    option = input('> ')
    if option == "1":
        print("Please provide the feature m/z, rt and CCS respectively (Ex. 150, 1, 150)")
        feat = input('> ')
        feat = [float(_) for _ in feat.split()]
        print("Please type the tolerance for m/z, rt, and CCS, respectively (Ex. '0.01 0.1 3.')\n\t"
              "* separated by spaces\n\t* CCS tolerance is a percentage, not an absolute value")
        tol = input('> ')
        tol = [float(_) for _ in tol.split()]
        mzs = [feat[0], tol[0]]
        rts = [feat[1], tol[1]]
        ccses = [feat[2], tol[2]]
        filtered = filter_d(mzs, rts, ccses, df)
        try:
            max_inten = max(filtered.iloc[0][3:])
        except Exception as e:
            print('! ERROR:', e)
            print("! ERROR: Unable find matching feature.")
            return False
        norm = []
        for i in range(3, len(filtered.iloc[0])):
            norm.append(max_inten / filtered.iloc[0][i])
        try:
            dset.normalize(np.asarray(norm))
            print('! INFO: Successfully normalized')
            return True
        except ValueError as ve:
            print('! ERROR:', ve)
            print("! ERROR: Unable to normalize. Please check the input and try again.")
            return False

    elif option == "2":
        print("Please provide a text file with the normalization values")
        path = input()
        norm = []
        try:
            with open(path) as fp:
                line = float(fp.readline())
                norm.append(line)
                while line:
                    try:
                        line = float(fp.readline())
                        norm.append(line)
                    except ValueError:
                        break
        except IOError as ioe:
            print('! ERROR:', ioe)
            print("! ERROR: Unable to normalize. Please check the path and try again.")
            return False
        try:
            dset.normalize(np.asarray(norm))
            print('! INFO: Successfully normalized')
            return True
        except ValueError as ve:
            print('! ERROR:', ve)
            print("! ERROR: Unable to normalize. Check the file and try again.")
            return False

    elif option == "back":
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


def calibrate_rt(dset):
    """
calibrate_rt
    description:
        Prompts the user with options to perform retention time calibration
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- finished calibrating retention time
"""
    def build_new_calibration():
        """ helper function that prompts user to enter RT calibrant information"""
        finished = False
        lipids, meas_rt, ref_rt = [], [], []
        while not finished:
            print("Adding lipid calibrants... Please enter a lipid to search for reference retention time")
            print("\t* example: 'PG(36:1)' *")
            print("\t'done' - finished adding calibrants")
            print("\t'cancel' - discard calibration")
            lipid = input('> ')
            if lipid == 'done':
                # check for lipid calibrants (at least one must be set)
                if len(lipids) < 1 or len(meas_rt) < 1 or len(ref_rt) < 1:
                    print('! INFO: did not create retention time calibration')
                    return None
                # build the actual RTCalibration and return it
                finished = True
            elif lipid == 'cancel':
                print('! INFO: discarding calibrants and cancelling retention time calibration creation')
                return None
            else:
                try:
                    parsed = parse_lipid(lipid)
                    fam = parsed['fa_mod'] if 'fa_mod' in parsed else None
                    rt_strict = get_ref_rt(parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat'], fa_mod=fam)
                    rt_nonstrict = get_ref_rt(parsed['lipid_class'], parsed['n_carbon'], parsed['n_unsat'],
                                              fa_mod=fam, strict=False)
                except:
                    rt_strict = None
                    rt_nonstrict = None
                if rt_strict is None and rt_nonstrict is None:
                    print('! ERROR: unable to find reference RT for lipid: {}'.format(lipid))
                else:
                    if rt_strict is not None:
                        print('! INFO: found reference RT for lipid: {} → {:.2f} min'.format(lipid, rt_strict))
                        rtr = rt_strict
                    else:
                        print('! WARNING: had to ignore unsaturations, RT matched only on lipid class and' \
                              ' fatty acid carbons')
                        print('! INFO: found reference RT for lipid: {} → {:.2f} min'.format(lipid, rt_nonstrict))
                        rtr = rt_nonstrict
                    # get the measured retention time for the lipid
                    print('Please enter the measured retention time for this lipid...')
                    rtm = float(input('> '))
                    lipids.append(lipid)
                    meas_rt.append(rtm)
                    ref_rt.append(rtr)
                    print('! INFO: added lipid calibrant: '
                          '{} measured RT → {:.2f} min reference RT → {:.2f} min'.format(lipid, rtm, rtr))
        rtc = RTCalibration(lipids, meas_rt, ref_rt)
        if rtc is not None:
            print('! INFO: successfully created a retention time calibration')
        else:
            print('! ERROR: failed to create a retention time calibration')
        return rtc

    print("Retention Time Calibration... Please choose an option")
    print("\t1. create new retention time calibration")
    print("\t2. view retention time calibration")
    print("\t3. clear current retention time calibration")
    print("\t'back' to go back")
    option = input('> ')

    if option == '1':
        # create a new RT calibration
        try:
            dset.rt_calibration = build_new_calibration()
        except Exception as e:
            print('! ERROR:', e)
            print('! ERROR: failed to create a new RT calibration.')
        # prompt again
        return False

    elif option == '2':
        if dset.rt_calibration is not None:
            # print the RT calibration information
            print(dset.rt_calibration)
        else:
            print('! INFO: there is no retention time calibration available')
        # prompt again
        return False

    elif option == '3':
        if dset.rt_calibration is None:
            print('! ERROR: there is no retention time calibration to clear')
            return False
        else:
            dset.rt_calibration = None
            print('! INFO: cleared existing retention time calibration')
            return False

    elif option == 'back':
        # go back
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        # prompt again
        return False


def batch_feature_selection(dset):
    """
batch_feature_selection
    description:
        Prompts the user to
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- success exporting features
"""
    print("Batch Feature Selection... Enter the name of the file with features to select")
    print("\texample: 'select_these_features.csv'")
    print('\t"back" to go back')
    inpath = input('> ')
    if inpath == "back":
        return True
    if not os.path.isfile(inpath):
        print('! ERROR: unable to find input file: "{}", check path and try again'.format(inpath))
        return False
    print("Enter the name of file to save the selected data in\n\texample: 'selected_features_raw.csv'")
    outpath = input('> ')
    print("Please type the search tolerance for m/z, rt, and CCS, respectively (Ex. '0.01 0.1 3.')\n\t"
          "* separated by spaces\n\t* CCS tolerance is a percentage, not an absolute value")
    tol = input('> ')
    tol = [float(_) for _ in tol.split()]
    success = False
    try:
        success = dset.select_feature_data(inpath, outpath, tolerance=tol)
        if not success:
            print('! ERROR: no matching features found in dataset')
            return False
        else:
            print('! INFO: selected features exported to: "{}"'.format(outpath))
    except Exception as e:
        print('! ERROR:', e)
        print('! ERROR: unable to select feature data ({})'.format(e))
    return True


def export(dset):
    """
export
    description:
        Prompts the user to save a excel file of the current lipidomics dataset:
            - The user specifies a path where to save the file and it creates an
              excel file with all the information that current lipidomics dataset
              has.
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
    returns:
        (bool) -- finished exporting data to excel file
"""
    print("Exporting data... Where would you like to save the file? \n\texample: 'jang_ho/results.xlsx'")
    print('\t"back" to go back')
    path = input('> ')
    if path == "back":
        return True
    try:
        dset.export_xlsx(path)
        print('! INFO: Successfully exported dataset to Excel spreadsheet: {}.'.format(path))
    except Exception as e:
        print('! ERROR:', e)
        print("! ERROR: Unable to export data.")
    return True


def main():
    """
main
    description:
        Main function for interactive interface. Prompts the user to choose actions ...
    parameters:
        () --
"""
    # keep retrying to load the Dataset until it returns something
    dset = load_dset()
    while dset is None:
        dset = load_dset()

    # instead of exiting, return None
    if dset == 'exit':
        return

    # create a pandas DataFrame
    label_df = DataFrame(dset.labels)
    int_df = DataFrame(dset.intensities)
    df = concat([label_df, int_df], axis=1, ignore_index=True, sort=False)

    # main execution loop
    while True:
        print("\nWhat would you like to do with this Dataset? ")
        print("\t1.  Manage Groups")
        print("\t2.  Filter Data")
        print("\t3.  Manage Statistics")
        print("\t4.  Make Plots")
        print("\t5.  Lipid Identification")
        print("\t6.  Normalize Intensities")
        print("\t7.  Calibrate Retention Time")
        print("\t8.  Overview of Dataset")
        print("\t9.  Batch Feature Selection")
        print("\t10. Export Current Dataset to Spreadsheet")
        print("\t11. Save Current Dataset to File")
        print('\t"exit" to quit the interface')
        option = input('> ')
        # Manage Groups
        if option == '1':
            finished = manage_groups(dset)
            while not finished:
                finished = manage_groups(dset)
        # Filter Data
        elif option == '2':
            finished = filter_data(dset)
            while not finished:
                finished = filter_data(dset)
        # Manage Statistics
        elif option == '3':
            finished = manage_statistics(dset)
            while not finished:
                finished = manage_statistics(dset)
        # Make Plots
        elif option == '4':
            finished = make_plots(dset)
            while not finished:
                finished = make_plots(dset)
        # Lipid Identification
        elif option == '5':
            identify_lipids(dset)
        # Normalize Intensities
        elif option == '6':
            finished = normalize_data(dset, df)
            while not finished:
                finished = normalize_data(dset, df)
        # Calibrate Retention Time
        elif option == '7':
            finished = calibrate_rt(dset)
            while not finished:
                finished = calibrate_rt(dset)
        # Overview of Dataset
        elif option == '8':
            print(dset)
        # Batch Feature Selection
        elif option == '9':
            batch_feature_selection(dset)
        # Export Current Dataset to Spreadsheet
        elif option == '10':
            export(dset)
        # Save Current Dataset to File
        elif option == '11':
            print("Saving Current Dataset to File... Please enter the full path and file name to save the Dataset "
                  "under.\n\t* .pickle file\n\t* no spaces in path)\n\texample: 'jang_ho/191120_bacterial_pos.pickle'")
            pickle_path = input('> ')
            try:
                dset.save_bin(pickle_path)
                print('! INFO: Dataset saved to file: "{}"'.format(pickle_path))
            except Exception as e:
                print('! ERROR:', e)
                print('! ERROR: Unable to save Dataset to file: "{}"'.format(pickle_path))
        elif option == 'exit':
            # instead of exiting return None
            return
        else:
            print('! ERROR: Unrecognized option: "{}"'.format(option))


if __name__ == "__main__":
    # main not setup to accept parameters yet...
    #main(sys.argv[1:])
    main()


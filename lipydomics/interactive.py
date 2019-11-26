"""
    lipydomics/interactive.py
    Jang Ho Cho and Dylan Ross

    description:

"""


import pandas as pd
import csv
import numpy as np
import os

from lipydomics.data import Dataset
from lipydomics.stats import add_anova_p, add_pca3, add_plsda, add_2group_corr
from lipydomics.plotting import (
    barplot_feature_bygroup, scatter_pca3_projections_bygroup, scatter_plsda_projections_bygroup,
    splot_plsda_pcorr_bygroup
)
from lipydomics.identification import add_feature_ids


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
            except:
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
        except ValueError:
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


def filter_d(mzs, rts, ccss, data):
    """ a helper function for filtering data
        given M/Z, RT, CCS ranges and a DataFrame containing data,
        find and returns all data within that range.
    """
    filtered = data[(data[0] < int(mzs[0]) + int(mzs[1])) & (data[0] > int(mzs[0]) - int(mzs[1])) &
                    (data[1] < int(rts[0]) + int(rts[1])) & (data[1] > int(rts[0]) - int(rts[1])) &
                    (data[2] < int(ccss[0]) + int(ccss[1])) & (data[2] > int(ccss[0]) - int(ccss[1]))]
    return filtered


""" To Do: S plot filter, """
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
    print('\t"back" to go back')
    option = input('> ')
    if option == "back":
        return True
    label_dat = dset.labels
    label_df = pd.DataFrame(label_dat)
    if option == "1":
        print("Please Provide M/Z range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
        mz = input('> ')
        print("Please Provide Retention Time range? (Ex. '1 1' <--- This would be 1 plus or minus 1)")
        rt = input('> ')
        print("Please Provide CCS range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
        ccs = input('> ')
        print("Which group would you like to choose? ('All' to select the whole data)")
        group = input('> ')
        try:
            if group == "All":
                cur_data = pd.DataFrame(dset.intensities)
            else:
                cur_data = dset.get_data_bygroup(group)
            int_df = pd.DataFrame(cur_data)
            cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
            mzs = mz.split()
            rts = rt.split()
            ccss = ccs.split()
            filtered = filter_d(mzs, rts, ccss, cur_df)
        except ValueError:
            return False
            print("! ERROR: Failed to filter data, please check your groups and try again")

    elif option == "2":
        print("Please provide the path of the file with batch-query information")
        path = input('> ')
        query = pd.read_csv(path)
        print("Which group would you like to choose? ('All' to select the whole data)")
        group = input('> ')
        try:
            if group == "All":
                cur_data = pd.DataFrame(dset.intensities)
                cur_df = pd.concat([label_df, cur_data], axis=1, ignore_index=True, sort=False)
            else:
                cur_data = dset.get_data_bygroup(group)
                int_df = pd.DataFrame(cur_data)
                cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
        except ValueError:
            print("! ERROR: Failed to filter data, please check your groups and try again")
            return False
        for index, row in query.iterrows():
            if index == 0:
                filtered = filter_d([row["m/z"], row["m/z_tol"]], [row["rt"], row["rt_tol"]],
                                    [row["ccs"], row["ccs_tol"]], cur_df)
            else:
                filtered = pd.concat([filtered,
                                      filter_d([int(row["m/z"]), int(row["m/z_tol"])],
                                               [row["rt"], row["rt_tol"]], [row["ccs"], row["ccs_tol"]],
                                               cur_df)])
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
        except:
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
        print('\t"back" to go back')
        opt2 = input('> ')
        # map options to stats functions
        stats_f_map = {'1': add_anova_p, '2': add_pca3, '3': add_plsda, '4': add_2group_corr}
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
                stats_f_map[opt2](dset, groups, normed=norm)
                print('! INFO: Applied new statistical analysis using groups: {}'.format(groups))
            except ValueError:
                print('! ERROR: Unable to perform statistical analysis, check group names and try again')
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
            #print(e)
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
    print("\t1. Bar plot feature(s) by group")
    print("\t2. Scatter PCA3 Projections by group")
    print("\t3. Scatter PLS-DA Projections by group")
    print("\t4. S-Plot PLSA-DA and Pearson correlation by group")
    print('\t"back" to go back')
    option = input('> ')

    # map the options to plotting functions
    plots_f_map = {'1': barplot_feature_bygroup, '2': scatter_pca3_projections_bygroup,
                   '3': scatter_plsda_projections_bygroup, '4': splot_plsda_pcorr_bygroup}
    if option in plots_f_map:
        print("Where would you like to save the plot? (default = current directory)")
        plot_dir = input('> ')
        print("Which groups would you like to plot?")
        groups = input('> ').split()
        print('Would you like to use normalized data? (y/N)')
        norm = input('> ') == 'y'
        if norm:
            # make sure there is actually normalized data in the Dataset
            if dset.normed_intensities is None:
                print('! ERROR: No normalization has been performed yet')
                return False
        if option == '1':
            try:
                # get the m/z, rt, and CCS of the feature
                print("Please enter the m/z, retention time, and CCS of the feature\n\t* separated by spaces\n\t*"
                      " example: '345.6789 1.23 123.45'")
                feature = input('> ').split()
                feat = [float(f) for f in feature]
                # get the m/z, rt, and CCS tolerance
                print("Please enter the search tolerances for m/z, retention time, and CCS\n\t* separated by spaces"
                      "\n\t* example: '0.05 0.5 1.0'")
                tolerance = input('> ').split()
                tol = [float(t) for t in tolerance]
                if plots_f_map[option](dset, groups, plot_dir, feat, normed=norm, tolerance=tol):
                    print('! INFO: Generated plot for groups: {}'.format(groups))
                else:
                    print('! ERROR: Could not match feature')
            except ValueError as ve:
                print(ve)
                print('! ERROR: Unable to generate plot using groups: {}'.format(groups))
        else:
            try:
                plots_f_map[option](dset, groups, plot_dir, normed=norm)
                print('! INFO: Generated plot for groups: {}'.format(groups))
            except KeyError as ke:
                print(ke)
                print('! ERROR: Required statistic(s) have not been computed for groups: {}'.format(groups))
            except ValueError as ve:
                print(ve)
                print('! ERROR: Unable to generate plot using groups: {}'.format(groups))
        return False

    elif option == 'back':
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


def identify_lipids(dset):
    """
identify lipids
    description:
        Prompts the user with options to perform lipid identification at a variety of levels of confidence
    parameters:
        dset (lipydomics.data.Dataset) -- lipidomics dataset instance
"""
    print("Identifying Lipids... Please enter the tolerances for m/z, retention time and CCS matching"
          "\n\t* separated by spaces\n\t* example: '0.05 0.5 1.0'")
    tolerance = input('> ').split()
    tol = [float(t) for t in tolerance]
    print("Please specify an identification level")
    print("\t'theo_mz' - match on theoretical m/z")
    print("\t'theo_mz_ccs' - match on theoretical m/z and CCS")
    print("\t'meas_mz_ccs' - match on measured m/z and CCS")
    print("\t'meas_mz_rt_ccs' - match on measured m/z, retention time, and CCS")
    print("\t'any' - try all criteria (highest confidence first)")
    print("\t'back' to go back")
    option = input('> ')

    if option in ['theo_mz', 'theo_mz_ccs', 'meas_mz_ccs', 'meas_mz_rt_ccs', 'any']:
        # make the identifications
        try:
            add_feature_ids(dset, tol, level=option)
            n_identified = len([_ for _ in dset.feat_ids if type(_) == list])
            print("! INFO: Lipid identifications added successfully ({} lipids identified)".format(n_identified))
        except ValueError:
            print("! ERROR: Unable to make lipid identifications")

    elif option == 'back':
        # just return None to go back
        return

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        # just return None to go back
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
        feat = input()
        feat = feat.split()
        print("Please type the tolerance for m/z, rt, and CCS, respectively (Ex. '0.1 0.1 0.01')")
        tol = input()
        tol = tol.split()
        mzs = [int(feat[0]), int(tol[0])]
        rts = [int(feat[1]), int(tol[1])]
        ccses = [int(feat[2]), int(tol[2])]
        filtered = filter_d(mzs, rts, ccses, df)
        try:
            max_inten = max(filtered.iloc[0][3:])
        except:
            print("! ERROR: Unable find matching feature.")
            return False
        norm = []
        for i in range(3, len(filtered.iloc[0])):
            norm.append(filtered.iloc[0][i] / max_inten)
        try:
            dset.normalize(np.asarray(norm))
            print('! INFO: Successfully normalized')
            return True
        except ValueError:
            print("! ERROR: Unable to normalize. Please check the constraints and try again.")
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
        except IOError:
            print("! ERROR: Unable to normalize. Please check the path and try again.")
            return False
        try:
            dset.normalize(np.asarray(norm))
            print('! INFO: Successfully normalized')
            return True
        except ValueError:
            print("! ERROR: Unable to normalize. Check the file and try again.")
            return False
    elif option == "back":
        return True
    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


def export(dset, df):
    """
    export
        description:
            Prompts the user to save a excel file of the current lipidomics dataset:
                - The user specifies a path where to save the file and it creates an
                  excel file with all the information that current lipidomics dataset
                  has.
        parameters:
            dset (lipydomics.data.Dataset) -- lipidomics dataset instance
            df (Pandas DataFrame) -- DataFrame version of lipidomics dataset
        returns:
            (bool) -- finished exporting data to excel file
    """
    print("Exporting data... Where would you like to save the file? \n\texample: 'jang_ho/results.xlsx'")
    print('\t"back" to go back')
    path = input('> ')
    if path == "back":
        return True
    writer = pd.ExcelWriter(path, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Data')
    for key in dset.stats:
        stats_df = pd.DataFrame(dset.stats[key])
        if "PCA3" in key and "loadings" in key:
            stats_df = stats_df.transpose()
        stats_df.to_excel(writer, sheet_name=key)
    m = 0
    feat_dict = {}
    label_df = pd.DataFrame(dset.labels)
    if dset.feat_ids:
        for feat in dset.feat_ids:
            if type(feat) is list:
                m = max(m, len(feat))

    for i in range(0, m):
        s = []
        for feat in dset.feat_ids:
            if type(feat) is list:
                try:
                    s.append(feat[i])
                except:
                    s.append("")
            else:
                if i == 0:
                    s.append(feat)
                else:
                    s.append("")
        feat_dict[i] = s
    if dset.normed_intensities is not None:
        norm_df = pd.DataFrame(dset.normed_intensities)
        norm_df = pd.concat([label_df, norm_df], axis=1, ignore_index=True, sort=False)

        norm_df.to_excel(writer, sheet_name="Normalized Intensities")
    if feat_dict:
        iden_df = pd.DataFrame(feat_dict)
        level_df = pd.DataFrame(dset.feat_id_levels)
        identification_df = pd.concat([label_df, level_df, iden_df], axis=1, ignore_index=True, sort=False)
        identification_df.to_excel(writer, sheet_name='Identification')
    try:
        writer.save()
        print('! INFO: Successfully exported excel version of the data to {}.'.format(path))
    except:
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

    # why?
    # create a pandas DataFrame
    label_df = pd.DataFrame(dset.labels)
    int_df = pd.DataFrame(dset.intensities)
    df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)

    # main execution loop
    while True:
        print("\nWhat would you like to do with this Dataset? ")
        print("\t1. Manage Groups")
        print("\t2. Filter Data")
        print("\t3. Manage Statistics")
        print("\t4. Make Plots")
        print("\t5. Lipid Identification")
        print("\t6. Normalize Intensities")
        print("\t7. Overview of Dataset")
        print("\t8. Export Current Dataset to Spreadsheet")
        print("\t9. Save Current Dataset to File")
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
        # Overview of Dataset
        elif option == '7':
            print(dset)
        # Export Current Dataset to Spreadsheet
        elif option == '8':
            export(dset, df)
        # Save Current Dataset to File
        elif option == '9':
            print("Saving Current Dataset to File... Please enter the full path and file name to save the Dataset "
                  "under.\n\t* .pickle file\n\t* no spaces in path)\n\texample: 'jang_ho/191120_bacterial_pos.pickle'")
            pickle_path = input('> ')
            try:
                dset.save_bin(pickle_path)
                print('! INFO: Dataset saved to file: "{}"'.format(pickle_path))
            except:
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


"""
    lipydomics/interactive.py
    Jang Ho Cho

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
    returns:
        (lipydomics.data.Dataset or None) -- lipidomics dataset instance
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
        # load the Dataset from .csv file
        dset = Dataset(csv_fname)
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
        exit()

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
              "\n\t- group name should not contain spaces\n\t- indices start at 0\n\t- example: 'A 1 3'")
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


def filter_data(dset):
    """
filter_data
    description:
        Prompts the user with options for filtering the data:
            -
            -
    parameters:
        () --
    returns:
        (bool) -- finished filtering data
"""
    def filter_d(mzs, rts, ccss, data):
        """ a helper function for filtering data """
        filtered = data[(data[0] < int(mzs[0]) + int(mzs[1])) & (data[0] > int(mzs[0]) - int(mzs[1])) &
                        (data[1] < int(rts[0]) + int(rts[1])) & (data[1] > int(rts[0]) - int(rts[1])) &
                        (data[2] < int(ccss[0]) + int(ccss[1])) & (data[2] > int(ccss[0]) - int(ccss[1]))]
        return filtered

    print('Filtering data... What would you like to do?')
    print("\t1. Single query")
    print("\t2. Batch query")
    print('\t"back" to go back')
    option = input('> ')

    return True

    if option == "1":
        pass
    label_dat = data.labels
    label_df = pd.DataFrame(label_dat)
    if option == "1":
        print(">> M/Z range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
        mz = input()
        print(">> Retention Time range? (Ex. '1 1' <--- This would be 1 plus or minus 1)")
        rt = input()
        print(">> CCS range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
        ccs = input()
        print(">> Which group would you like to choose? ('All' to select the whole data)")
        group = input()
        if group == "b":
            pass
        try:
            if group == "All":
                cur_data = pd.DataFrame(data.intensities)
            else:
                cur_data = data.get_data_bygroup(group)
            int_df = pd.DataFrame(cur_data)
            cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
            mzs = mz.split()
            rts = rt.split()
            ccss = ccs.split()
            filtered = filter_d(mzs, rts, ccss, cur_df)
            print(filtered)
        except ValueError:
            print(" >> That group has not been assigned")

    elif option == "2":
        print(">> Which group would you like to choose? ('All' to select the whole data)")
        group = input()
        if group == "b":
            pass
        try:
            if group == "All":
                cur_data = pd.DataFrame(data.intensities)
            else:
                cur_data = data.get_data_bygroup(group)
                int_df = pd.DataFrame(cur_data)
                cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
        except ValueError:
            print(">> That group has not been assigned")
        print(">> Path of the file with batch-query information")
        path = input()
        query = pd.read_csv(path)
        for index, row in query.iterrows():
            if index == 0:
                filtered = filter_d([row["m/z"], row["m/z_tol"]], [row["rt"], row["rt_tol"]],
                                    [row["ccs"], row["ccs_tol"]], cur_df)
            else:
                filtered = pd.concat([filtered,
                                      filter_d([int(row["m/z"]), int(row["m/z_tol"])],
                                               [row["rt"], row["rt_tol"]], [row["ccs"], row["ccs_tol"]],
                                               cur_df)])
    print(filtered)
    print(">> Data filter success. Would you like to download the result as csv? (y/n)")
    option = input()
    if option == "y":
        print(">> Please specify the path you want to save the csv file")
        path = input()
        export_csv = filtered.to_csv('results.csv',
                                     index=None, header=False)

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
    print("\t3. Download CSV file of computed Statistics")
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
            print('Please enter group names for this analysis, separated by spaces')
            groups = input('> ').split()
            try:
                stats_f_map[opt2](dset, groups, normed=norm)
            except ValueError:
                print('! ERROR: Unable to perform statistical analysis, check group names and try again')
        elif opt2 == 'back':
            pass
        else:
            print('! ERROR: unrecognized option: "{}"'.format(opt2))
        return False

    elif option == '2':
        return False

    elif option == '3':
        return False

    elif option == 'back':
        return True

    else:
        print('! ERROR: unrecognized option: "{}"'.format(option))
        return False


"""if option == "1":
                print(">> What would you like to do?")
                print("1. Anova-P")
                print("2. PCA3")
                print("3. PLS-DA")
                print("4. Two Group Correlation")
                option = input()
                if option == "b":
                    continue
                print(">> Which groups would you like to perform chosen Statistic on?")
                group = input()
                print(">> Would you like to use normalized data? (y/n)")
                norm = input()
                if norm == "y":
                    norm = True
                else:
                    norm = False
                if group == "b":
                    continue
                group = group.split()
                if option == "1":
                    try:
                        add_anova_p(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "2":
                    try:aaaa
                        add_pca3(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "3":
                    try:
                        add_plsda(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "4":
                    try:
                        add_2group_corr(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
            elif option == "2":
                print(data.stats)
            elif option == "3":
                path = ""
                df = pd.DataFrame.from_dict(data.stats)
                export = df.to_csv(r'C:/Users/narsi/Desktop/results.csv')"""


def main():
    """
main
    description:
        Main function for interactive interface. Prompts the user to choose actions ...
    parameters:
        () --
"""
    # keep retrying to load the Dataset until it works (or the exit command happens)
    dset = load_dset()
    while dset is None:
        dset = load_dset()

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
        print("\t8. Download current data as CSV")
        print("\t9. Save current Dataset")
        print('\t"exit" to quit the interface')
        option = input('> ')
        # Manage groups
        if option == '1':
            finished = manage_groups(dset)
            while not finished:
                finished = manage_groups(dset)
        elif option == '2':
            finished = filter_data(dset)
            while not finished:
                finished = filter_data(dset)
        elif option == '3':
            finished = manage_statistics(dset)
            while not finished:
                finished = manage_statistics(dset)
        elif option == '4':
            pass
        elif option == '5':
            pass
        elif option == '6':
            pass
        elif option == '7':
            print(dset)
        elif option == '8':
            pass
        elif option == '9':
            pass
        elif option == 'exit':
            exit()
        else:
            print('! ERROR: Unrecognized option: "{}"'.format(option))


"""
        

        elif option == "2":
            print("1. Single query")
            print("2. Batch query")
            option = input()
            if option == "b":
                continue
            label_dat = data.labels
            label_df = pd.DataFrame(label_dat)
            if option == "1":
                print(">> M/Z range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
                mz = input()
                print(">> Retention Time range? (Ex. '1 1' <--- This would be 1 plus or minus 1)")
                rt = input()
                print(">> CCS range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
                ccs = input()
                print(">> Which group would you like to choose? ('All' to select the whole data)")
                group = input()
                if group == "b":
                    continue
                try:
                    if group == "All":
                        cur_data = pd.DataFrame(data.intensities)
                    else:
                        cur_data = data.get_data_bygroup(group)
                    int_df = pd.DataFrame(cur_data)
                    cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
                    mzs = mz.split()
                    rts = rt.split()
                    ccss = ccs.split()
                    filtered = filter_d(mzs, rts, ccss, cur_df)
                    print(filtered)
                except ValueError:
                    print(" >> That group has not been assigned")

            elif option == "2":
                print(">> Which group would you like to choose? ('All' to select the whole data)")
                group = input()
                if group == "b":
                    continue
                try:
                    if group == "All":
                        cur_data = pd.DataFrame(data.intensities)
                    else:
                        cur_data = data.get_data_bygroup(group)
                        int_df = pd.DataFrame(cur_data)
                        cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
                except ValueError:
                    print(">> That group has not been assigned")
                print(">> Path of the file with batch-query information")
                path = input()
                query = pd.read_csv(path)
                for index, row in query.iterrows():
                    if index == 0:
                        filtered = filter_d([row["m/z"], row["m/z_tol"]], [row["rt"], row["rt_tol"]],
                                            [row["ccs"], row["ccs_tol"]], cur_df)
                    else:
                        filtered = pd.concat([filtered,
                                              filter_d([int(row["m/z"]), int(row["m/z_tol"])],
                                                       [row["rt"], row["rt_tol"]], [row["ccs"], row["ccs_tol"]],
                                                       cur_df)])
            print(filtered)
            print(">> Data filter success. Would you like to download the result as csv? (y/n)")
            option = input()
            if option == "y":
                print(">> Please specify the path you want to save the csv file")
                path = input()
                export_csv = filtered.to_csv('results.csv',
                                             index=None, header=False)

        elif option == "3":
            print("1. Compute Statistics")
            print("2. View Statistics")
            print("3. Download CSV file of computed Statistics")
            option = input()
            if option == "1":
                print(">> What would you like to do?")
                print("1. Anova-P")
                print("2. PCA3")
                print("3. PLS-DA")
                print("4. Two Group Correlation")
                option = input()
                if option == "b":
                    continue
                print(">> Which groups would you like to perform chosen Statistic on?")
                group = input()
                print(">> Would you like to use normalized data? (y/n)")
                norm = input()
                if norm == "y":
                    norm = True
                else:
                    norm = False
                if group == "b":
                    continue
                group = group.split()
                if option == "1":
                    try:
                        add_anova_p(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "2":
                    try:aaaa
                        add_pca3(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "3":
                    try:
                        add_plsda(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "4":
                    try:
                        add_2group_corr(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
            elif option == "2":
                print(data.stats)
            elif option == "3":
                path = ""
                df = pd.DataFrame.from_dict(data.stats)
                export = df.to_csv(r'C:/Users/narsi/Desktop/results.csv')
        elif option == "4":
            print("1. Bar plot feature by group")
            print("2. Scatter PCA3 Projections by group")
            print("3. Scatter PLS-DA Projections by group")
            print("4. S-Plot PLSA-DA and Pearson correlation by group")

            option = input()
            if option == "b":
                continue
            print(">> Where would you like to save the plot?")
            path = ""
            print(">> Which group would you like to plot?")
            group = input()
            if group == "b":
                continue
            group = group.split()
            print(">> Would you like to use normalized data? (y/n)")
            norm = input()
            if norm == "y":
                norm = True
            else:
                norm = False
            if option == "1":
                print(">> Feature Range? (Type M/Z RT CCS in this order)")
                feature = input()
                feature = list(map(float, feature.split()))
                barplot_feature_bygroup(data, group, path, norm, feature)
            if option == "2":
                try:
                    scatter_pca3_projections_bygroup(data, group, path, norm)
                except KeyError:
                    print("Statistic not yet computed")
            if option == "3":
                try:
                    scatter_plsda_projections_bygroup(data, group, path, norm)
                except KeyError:
                    print("Statistic not yet computed")
            if option == "4":
                try:
                    splot_plsda_pcorr_bygroup(data, group, path, norm)
                except KeyError:
                    print("Statistic not yet computed")
        elif option == "5":
            print(">> Please type the tolerance for m/z, rt, and CCS, respectively (Ex. '0.1 0.1 0.01')")
            feature = input()
            if feature == "b":
                continue
            feature = list(map(float, feature.split()))
            print(">> Level? (theo_mz', 'meas_mz_rt_ccs', 'any')")
            level = input()
            if level == "b":
                continue
            try:
                add_feature_ids(data, feature, level)
                print("Identification added successfully")
            except ValueError:
                print("Something went wrong!")

        elif option == "6":
            print("1. Internal")
            print("2. External")
            option = input()
            if option == "1":
                print(">> Please provide the feature m/z, rt and CCS respectively (Ex. 150, 1, 150)")
                feat = input()
                feat = feat.split()
                print(">> Please type the tolerance for m/z, rt, and CCS, respectively (Ex. '0.1 0.1 0.01')")
                tol = input()
                tol = tol.split()
                mzs = [int(feat[0]), int(tol[0])]
                rts = [int(feat[1]), int(tol[1])]
                ccses = [int(feat[2]), int(tol[2])]
                filtered = filter_d(mzs, rts, ccses, df)
                max_inten = max(filtered.iloc[0][3:])
                norm = []
                for i in range(3, len(filtered.iloc[0])):
                    norm.append(filtered.iloc[0][i] / max_inten)
                try:
                    data.normalize(np.asarray(norm))
                    print("Successfully normalized")
                except ValueError:
                    print("Something went wrong")

            elif option == "2":
                print("    >> Please provide a text file with the normalization values")
                path = input()
                norm = []
                with open(path) as fp:
                    line = float(fp.readline())
                    norm.append(line)
                    while line:
                        try:
                            line = float(fp.readline())
                            norm.append(line)
                        except ValueError:
                            break
                try:
                    data.normalize(np.asarray(norm))
                    print("Successfully normalized")
                except ValueError:
                    print("Something went wrong")

        elif option == "7":
            print(data)

        elif option == "8":
            print("Where would you like to save?")
            path = input()
            if path == "b":
                continue
            writer = pd.ExcelWriter('results.xlsx', engine='xlsxwriter')
            df.to_excel(writer, sheet_name='Data')
            for key in data.stats:
                stats_df = pd.DataFrame(data.stats[key])
                if "PCA3" in key and "loadings" in key:
                    stats_df = stats_df.transpose()
                stats_df.to_excel(writer, sheet_name=key)
            m = 1
            feat_dict = {}
            for feat in data.feat_ids:
                if type(feat) is list:
                    m = max(m, len(feat))
            for i in range(0, m):
                s = []
                for feat in data.feat_ids:
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
            if data.normed_intensities is not None:
                norm_df = pd.DataFrame(data.normed_intensities)
                norm_df.to_excel(writer, sheet_name="Normalized Intensities")
            iden_df = pd.DataFrame(feat_dict)
            level_df = pd.DataFrame(data.feat_id_levels)
            pd.concat([level_df, iden_df], axis=1, ignore_index=False, sort=False)
            identification_df = pd.concat([level_df, iden_df], axis=1, ignore_index=True, sort=False)
            identification_df.to_excel(writer, sheet_name='Identification')
            writer.save()
        elif option == "9":
            print("Where do you want to save your data set?")
            path = input()
            if path == "b":
                continue
            data.save_bin(path)
            print("File saved successfully")
        elif option == "exit":
            c = False
            
"""

if __name__ == "__main__":
    # main not setup to accept parameters yet...
    #main(sys.argv[1:])
    main()


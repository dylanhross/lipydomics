from data import Dataset
import stats
import plotting
import identification
import pandas as pd
import csv
import numpy as np
import sys


def main():
    print("1. Make a new Dataset")
    print("2. Load previous Dataset")
    option = input()
    c = True

    def filter_d(mzs, rts, ccss, data):
        filtered = data[(data[0] < int(mzs[0]) + int(mzs[1])) & (data[0] > int(mzs[0]) - int(mzs[1])) &
                        (data[1] < int(rts[0]) + int(rts[1])) & (data[1] > int(rts[0]) - int(rts[1])) &
                        (data[2] < int(ccss[0]) + int(ccss[1])) & (data[2] > int(ccss[0]) - int(ccss[1]))]
        return filtered

    if option == "1":
        print("Please enter the path of the csv file you want to work with: ")
        file = input()
        try:
            print("Does this file have headers? y/n")
            ans = input()
            if ans == "y":
                print("Would you like to auto-make groups? (y/n)")
                ans = input()
                if ans == "y":
                    with open(file, newline='') as f:
                        reader = csv.reader(f)
                        header = next(reader)
                    header = header[3:]
                    group_map = {}
                    for i in range(len(header)):
                        if header[i] not in group_map:
                            group_map[header[i]] = [i]
                        else:
                            group_map[header[i]] = group_map[header[i]] + [i]
                    data = Dataset(file)
                    data.assign_groups(group_map)
                else:
                    data = Dataset(file)
            else:
                data = Dataset(file)
            print("Data loaded successfully")
        except IOError:
            print(">> Error. Make sure the path specified is correct and the file exists")
            c = False
    elif option == "2":
        print("Please provide the path of data file you want to load")
        path = input()
        try:
            data = Dataset.load_bin(path)
            print("Data loaded successfully")
        except IOError:
            print(">> Error. Make sure the path specified is correct and the file exists")
            c = False
    if c:
        print(data)
        label_df = pd.DataFrame(data.labels)
        int_df = pd.DataFrame(data.intensities)
        df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
    while c:
        print("What would you like to do with the data? ")
        print("1. Manage groups")
        print("2. Filter Data")
        print("3. Manage Statistics")
        print("4. Make Plots")
        print("5. Identification")
        print("6. Normalize")
        print("7. Overview of Dataset")
        print("8. Download current data as CSV")
        print("9. Save current Dataset")
        print("'exit' to quit the interface")
        option = input()
        if option == "1":
            print(" 1. Assign Groups")
            print(" 2. View Assigned Groups")
            print(" 3. Get Data By Group(s)")
            option = input()
            if option == "b":
                continue
            if option == "1":
                cont = True
                while cont:
                    print(">> Please provide a name for a group and its indices in order of Name > Starting index > "
                          "Ending index.\nExample: 'A 1 3'")
                    group = input()
                    if group == "b":
                        break
                    group = group.split()
                    d = {group[0]: range(int(group[1]), int(group[2]) + 1)}
                    try:
                        data.assign_groups(d)
                        print("     >> Group assigned successfully")
                        print("     >> Would you like to assign more groups? (y/n)")
                        a = input()
                        if a == "n":
                            cont = False
                    except ValueError:
                        print(
                            ">> Error. Make sure the provided indices are correct and the format of the input is"
                            " like the provided example")

            if option == "2":
                print(" " + data.group_indices)

            if option == "3":
                print(">> Which group would you like to view?")
                name = input()
                print(data.get_data_bygroup(name))

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
                        stats.add_anova_p(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "2":
                    try:
                        stats.add_pca3(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "3":
                    try:
                        stats.add_plsda(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
                if option == "4":
                    try:
                        stats.add_2group_corr(data, group, norm)
                        print(">> Statistic added successfully")
                    except ValueError:
                        print(">> Something went wrong")
            elif option == "2":
                print(data.stats)
            elif option == "3":
                path = ""
                df = pd.DataFrame.from_dict(data.stats)
                export = df.to_csv(r'C:\Users\narsi\Desktop\results.csv')
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
                plotting.barplot_feature_bygroup(data, group, path, norm, feature)
            if option == "2":
                try:
                    plotting.scatter_pca3_projections_bygroup(data, group, path, norm)
                except KeyError:
                    print("Statistic not yet computed")
            if option == "3":
                try:
                    plotting.scatter_plsda_projections_bygroup(data, group, path, norm)
                except KeyError:
                    print("Statistic not yet computed")
            if option == "4":
                try:
                    plotting.splot_plsda_pcorr_bygroup(data, group, path, norm)
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
                identification.add_feature_ids(data, feature, level)
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


if __name__ == "__main__":
    main(sys.argv[1:])


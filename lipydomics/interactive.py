import data as Dataset
import stats
import plotting
import pandas as pd

if __name__ == "__main__":
    print("Please enter the path of the csv file you want to work with: ")
    file = input()
    try:
        data = Dataset(file)
        print("Data loaded successfully")
        print(data)
        c = True
        while c:
            print("What would you like to do with the data? ")
            print("1. Manage groups")
            print("2. Filter Data")
            print("3. Manage Statistics")
            print("4. Make Plots")
            print("'exit' to quit the interface")
            option = input()
            if option == "1":
                print(" 1. Assign Groups")
                print(" 2. View Assigned Groups")
                print(" 3. Get Data By Group(s)")
                option = input()
                if option == "1":
                    cont = True
                    while cont:
                        print("     >> Please provide a name for a group and its indices in order of Name > Starting index > "
                              "Ending index.\nExample: 'A 1 3'")
                        group = input()
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
                                "       >> Error. Make sure the provided indices are correct and the format of the input is"
                                " like the provided example")

                if option == "2":
                    print(" " + data.group_indices)

                if option == "3":
                    print(" >> Which group would you like to view?")
                    name = input()
                    print(data.get_data_bygroup(name))

            elif option == "2":
                print(" 1. Single query")
                print(" 2. Batch query")
                option = input()
                label_dat = data.labels
                label_df = pd.DataFrame(label_dat)
                def filter_d(mz, rt, ccs, data):
                    filtered = data[(cur_df[0] < int(mzs[0]) + int(mzs[1])) & (cur_df[0] > int(mzs[0]) - int(mzs[1])) &
                                      (cur_df[1] < int(rts[0]) + int(rts[1])) & (cur_df[1] > int(rts[0]) - int(rts[1])) &
                                      (cur_df[2] < int(ccss[0]) + int(ccss[1])) & (cur_df[2] > int(ccss[0]) - int(ccss[1]))]
                    return filtered
                if option == "1":
                    print(" >> M/Z range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
                    mz = input()
                    print(" >> Retention Time range? (Ex. '1 1' <--- This would be 1 plus or minus 1)")
                    rt = input()
                    print(" >> CCS range? (Ex. '150 10'  <--- This would be 150 plus or minus 10)")
                    ccs = input()
                    print(" >> Which group would you like to choose? ('All' to select the whole data)")
                    group = input()
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
                    print(" >> Which group would you like to choose? ('All' to select the whole data)")
                    group = input()
                    try:
                        if group == "All":
                                cur_data = pd.DataFrame(data.intensities)
                        else:
                            cur_data = data.get_data_bygroup(group)
                            int_df = pd.DataFrame(cur_data)
                            cur_df = pd.concat([label_df, int_df], axis=1, ignore_index=True, sort=False)
                    except ValueError:
                        print(" >> That group has not been assigned")
                    print(" >> Path of the file with batch-query information")
                    path = input()
                    query = pd.read_csv(path)
                    for index, row in query.iterrows():
                        if index == 0:
                            filtered = filter_d([row["m/z"], row["m/z_tol"]], [row["rt"], row["rt_tol"]], [row["ccs"], row["ccs_tol"]], cur_df)
                        else:
                            filtered = pd.concat([filtered,
                                                  filter_d([int(row["m/z"]), int(row["m/z_tol"])], [row["rt"], row["rt_tol"]], [row["ccs"], row["ccs_tol"]], cur_df)])
                print(filtered)
                print(" >> Data filter success. Would you like to download the result as csv? (y/n)")
                option = input()
                if option == "y":
                    print("    >> Please specify the path you want to save the csv file")
                    path = input()
                    export_csv = filtered.to_csv ('results.csv',
                                            index = None, header=False)

            elif option == "3":
                print(" 1. Compute Statistics")
                print(" 2. View Statistics")
                print(" 3. Download CSV file of computed Statistics")
                option = input()
                if option == "1":
                    print("     >> What would you like to do?")
                    print("        1. Anova-P")
                    print("        2. PCA3")
                    print("        3. PLS-DA")
                    print("        4. Two Group Correlation")
                    option = input()
                    print("         >> Which groups would you like to perform chosen Statistic on?")
                    group = input()
                    group = group.split()
                    if option == "1":
                        try:
                            stats.add_anova_p(data, group)
                            print("            >> Statistic added successfully")
                        except ValueError:
                            print("            >> Something went wrong")
                    if option == "2":
                        try:
                            stats.add_pca3(data, group)
                            print("            >> Statistic added successfully")
                        except ValueError:
                            print("            >> Something went wrong")
                    if option == "3":
                        try:
                            stats.add_plsda(data, group)
                            print("            >> Statistic added successfully")
                        except ValueError:
                            print("            >> Something went wrong")
                    if option == "4":
                        try:
                            stats.add_2group_corr(data, group)
                            print("            >> Statistic added successfully")
                        except ValueError:
                            print("            >> Something went wrong")
                elif option == "2":
                    print(data.stats)
                elif option == "3":
                    path = ""
                    df = pd.DataFrame.from_dict(data.stats)
                    export = df.to_csv(r'C:\Users\narsi\Desktop\results.csv')
            elif option == "4":
                print(" 1. Bar plot feature by group")
                print(" 2. Scatter PCA3 Projections by group")
                print(" 3. Scatter PLS-DA Projections by group")
                print(" 4. S-Plot PLSA-DA and Pearson correlation by group")

                option = input()
                print("    >> Where would you like to save the plot?")
                path = ""
                print("    >> Which group would you like to plot?")
                group = input()
                group = group.split()
                if option == "1":
                    print("    >> Feature Range? (Type M/Z RT CCS in this order)")
                    feature = input()
                    feature = list(map(float, feature.split()))
                    plotting.barplot_feature_bygroup(data, group, path, feature)
                if option == "2":
                    try:
                        plotting.scatter_pca3_projections_bygroup(data, group, path)
                    except KeyError:
                        print("Statistic not yet computed")
                if option == "3":
                    try:
                        plotting.scatter_plsda_projections_bygroup(data, group, path)
                    except KeyError:
                        print("Statistic not yet computed")
                if option == "4":
                    try:
                        plotting.splot_plsda_pcorr_bygroup(data, group, path)
                    except KeyError:
                        print("Statistic not yet computed")

            elif option == "exit":
                c = False

    except IOError:
        print(">> Error. Make sure the path specified is correct and the file exists")
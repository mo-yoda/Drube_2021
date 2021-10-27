import pandas as pd
import os
import numpy as np


def main():
    path = "C:/path/to/folder/"
    path_laptop = ""

    export_csv(path_to_folder=path)
    # export_csv(path_to_folder=path_laptop)


def import_df(path_to_folder):
    csv_list = extract_files_to_list(path_to_folder, path_bool=True)
    csv_name_list = extract_files_to_list(path_to_folder, path_bool=False)

    # list of csv files which will be processed
    print(csv_list)

    # initialising the list of formatted tables
    # curve -> concentration response data
    # baseline -> baseline and stimulated BRET ratio
    dataframe_list_curve = []
    curve_name_list = []
    dataframe_list_baseline = []
    baseline_name_list = []

    for i, csv in enumerate(csv_list):
        data = pd.read_csv(csv, sep=",")
        # print(data)

        # if first condition is dQ + EV -> concentration response data
        if data.columns[1] == "Q-GRK + EV":
            dataframe_list_curve.append(data)
            curve_name_list.append(csv_name_list[i])
        # if first condition is "baseline" -> baseline and stimulated BRET data
        if data.columns[1] == "baseline":
            dataframe_list_baseline.append(data)
            baseline_name_list.append(csv_name_list[i])

        # print(baseline_name_list)
        # print(curve_name_list)

    return curve_name_list, dataframe_list_curve, baseline_name_list, dataframe_list_baseline

# formatting of concentration response data
def curve_format_df(path_to_folder):
    curve_name_list, dataframe_list_curve, baseline_name_list, dataframe_list_baseline = import_df(path_to_folder)

    # initialising the list of formatted tables
    curve_data_finished = []

    for data in dataframe_list_curve:
        # calculate n (3x technical replicates per n)
        n = int((len(data.columns) - 1) / 18)

        # delete first columns containing the applied ligand concentration
        data = data.drop(data.columns[0], axis=1)

        # delete all rows except the highest stimulation
        data = data.drop(data.index[1:len(data.index)])

        # exclude values (NaN) which were excluded in prism
        # excluded values are marked with *
        # to do this, transform pandas data frame to numpy array
        for i, value in enumerate(data.iloc[0]):
            try:
                data.iloc[:, i] = pd.to_numeric(value)
            except ValueError:
                data.iloc[:, i] = np.nan

        # calculate mean of technical triplicates for each n
        mean_data = pd.concat([data.iloc[:, i:i + 3].mean(axis=1) for i in range(0, len(data.columns), 3)], axis=1)
        mean_data = mean_data.transpose()

        # column names + condition
        mean_data = mean_data.rename(columns={0: "fold_change"})
        column = ["dQ", "GRK2", "GRK3", "GRK5", "GRK6", "CON"]
        condition = []

        for i in column:
            for j in range(n):
                condition += [i]

        mean_data["condition"] = condition

        # print(mean_data)

        # list of csv files which will be processed
        curve_data_finished.append(mean_data)
        print(curve_name_list)
    return curve_data_finished, curve_name_list


def baseline_format_df(path_to_folder):
    curve_name_list, dataframe_list_curve, baseline_name_list, dataframe_list_baseline = import_df(path_to_folder)

    # initialising the list of formatted tables
    baseline_data_finished = []

    for data in dataframe_list_baseline:
        # calculate n
        n = int((len(data.columns) - 1) / 2 / 24)

        # delete first column with condition names (contains umlauts)
        data = data.drop(data.columns[0], axis=1)

        # error if one row or column to many was imported (all NaN)
        # to solve this:
        data = data.dropna(axis=0, how="all")

        if len(data.columns) > n * 24 * 2 + 1:
            temp = len(data.columns) - n * 24 * 2
            for i in range(temp):
                data = data.drop(data.columns[len(data.columns) - 1], axis=1)

        data = data.rename(index={0: "dQ", 1: "GRK2", 2: "GRK3", 3: "GRK5", 4: "GRK6", 5: "CON"})

        # exclude values (NaN) which were excluded in prism
        # excluded values are marked with *
        # to do this, transform pandas data frame to numpy array
        data = data.to_numpy()

        for row in range(0, 6):
            for i, value in enumerate(data[row]):
                try:
                    data[row, i] = pd.to_numeric(value)
                except ValueError:
                    data[row, i] = np.nan

        # preparation of column and rownames for retransform numpy array to pandas data frame
        header = ["baseline", "stimulated"]
        column_names = []

        for i in header:
            for j in range(n * 24):
                column_names += [i]

        # retransform numpy array to pandas data frame
        data = pd.DataFrame(data=data, index=["dQ", "GRK2", "GRK3", "GRK5", "GRK6", "CON"], columns=column_names)

        # calculate mean for each n

        # mean of baselines (24x technical replicates per n)
        mean_baseline_data = pd.concat([data.iloc[:, i:i + 24].mean(axis=1) for i in range(0, n * 24, 24)], axis=1)
        # print(mean_baseline_data)

        # mean of stimulates data (3x technical replicates per n)
        mean_stimulated_data = pd.concat([data.iloc[:, i:i + 3].mean(axis=1) for i in range(n * 24, n * 24 + n * 3, 3)],
                                         axis=1)
        # print(mean_stimulated_data)

        # create column for id
        id = []
        for j in range(len(header)):
            for i in range(1, (len(data.index) * n) + 1):
                temp = str(i)
                id += [temp]

        # create column for condition
        condition = []
        for j in range(len(header)):
            for i in data.index:
                for r in range(n):
                    condition += [i]

        # create column for state (base or stim)
        state = []
        for j in header:
            for i in range(n * len(data.index)):
                state += [j]

        # create column for BRET ratio
        BRET = []

        # for loop extracts BRET ratios row wise and saves them in array
        # baseline BRET ratios
        for row in range(len(data.index)):
            for value in mean_baseline_data.iloc[row]:
                BRET += [value]

        # append stimulated BRET ratios
        for row in range(len(data.index)):
            for value in mean_stimulated_data.iloc[row]:
                BRET += [value]

        mean_data = pd.DataFrame(id)
        mean_data = mean_data.rename(columns={0: "id"})
        mean_data["condition"] = condition
        mean_data["state"] = state
        mean_data["BRET"] = BRET
        # print(mean_data)

        # find NaN by row
        is_NaN = mean_data.isnull()
        row_has_NaN = is_NaN.any(axis=1)
        rows_with_NaN = mean_data[row_has_NaN]
        # print(rows_with_NaN)

        # if rows containing NaNs were found, delete all row with this id
        # (if id = 1 in baseline is NaN, corresponding stimulated values are also deleted)
        if rows_with_NaN.empty == False:
            for id in rows_with_NaN["id"]:
                mean_data = mean_data[mean_data["id"] != id]

        # print(mean_data)

        # list of processed data frames
        baseline_data_finished.append(mean_data)
    return baseline_data_finished, baseline_name_list


def export_csv(path_to_folder):
    curve_data_finished, curve_name_list = curve_format_df(path_to_folder)
    print(curve_name_list)
    print(curve_data_finished)
    for i, table in enumerate(curve_data_finished):
        # print(table)
        curve_results_folder = path_to_folder + "Formatted curve data/"
        if not os.path.isdir(curve_results_folder):
            os.makedirs(curve_results_folder)
        table.to_csv(curve_results_folder + curve_name_list[i][0:len(curve_name_list[i]) - 4] + "_format.csv",
                     index=False)

    baseline_data_finished, baseline_name_list = baseline_format_df(path_to_folder)

    for i, table in enumerate(baseline_data_finished):
        baseline_results_folder = path_to_folder + "Formatted baseline data/"
        if not os.path.isdir(baseline_results_folder):
            os.makedirs(baseline_results_folder)
        table.to_csv(baseline_results_folder + baseline_name_list[i][0:len(baseline_name_list[i]) - 4] + "_format.csv",
                     index=False)

    return


def extract_files_to_list(path_to_folder, path_bool):
    """
    This functions checks a directory for all files of a certain dataype and returns a list of all found files
    :param path_to_folder: string
        Path to directory to be searched
    :param path_bool: boolean
        append path to filename if set to TRUE, only filename if FALSE
    :return: list
        Contains paths to all files or all filenames in searched directory
    """
    new_list = []
    for filename in os.listdir(path_to_folder):
        if filename.endswith(".csv"):
            if path_bool:
                new_list.append(os.path.join(path_to_folder, filename))
            else:
                new_list.append(filename)
        else:
            continue
    return new_list


if __name__ == '__main__':
    main()

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
    dataframe_list = []
    name_list = []

    for i, csv in enumerate(csv_list):
        data = pd.read_csv(csv, sep=",")
        # print(data)

        dataframe_list.append(data)
        name_list.append(csv_name_list[i])

    return name_list, dataframe_list

def format_df(path_to_folder):
    name_list, dataframe_list = import_df(path_to_folder)

    # initialising the list of processed tables
    data_finished = []
    t = 0
    for data in dataframe_list:
        print(name_list[t])
        t += 1
        # calculate n
        n = int((len(data.columns) - 1) / 2 / 24)

        # delete first column with condition names (contains umlauts)
        data = data.drop(data.columns[0], axis=1)

        data = data.rename(index={0: "barr", 1: "Goe", 2: "dFLR"})
        print(data)

        # exclude values which were excluded in prism and therefore in bargraph
        # excluded values are marked with *
        # to do this, transform pandas data frame to numpy array
        data = data.to_numpy()

        # for row in range(0, 6):
        #     for i, value in enumerate(data[row]):
        #         try:
        #             data[row, i] = pd.to_numeric(value)
        #         except ValueError:
        #             data[row, i] = np.nan

        # preparation of column and rownames for retransform numpy array to pandas data frame
        header = ["baseline", "stimulated"]
        column_names = []

        for i in header:
            if len(data[1]) < (2*256):
                for j in range(n * 24):
                    column_names += [i]
            else:
                for j in range(256):
                    column_names += [i]

        # retransform numpy array to pandas data frame
        data = pd.DataFrame(data=data, index=["barr", "Goe", "dFLR"], columns=column_names)
        print(data)

        # calculate mean of each triplicate

        if len(column_names) < (2*256):
            # mean of baselines (24x technical replicates per n)
            mean_baseline_data = pd.concat([data.iloc[:, i:i + 24].mean(axis=1) for i in range(0, n * 24, 24)], axis=1)

            # mean of stimulates data (3x technical replicates per n)
            mean_stimulated_data = pd.concat([data.iloc[:, i:i + 3].mean(axis=1) for i in range(n * 24, n * 24 + n * 3, 3)],
                                             axis=1)
        if len(column_names) == (2*256):
            mean_baseline_data = pd.concat([data.iloc[:, i:i + 24].mean(axis=1) for i in range(0, 9*24, 24)], axis=1)
            mean_baseline_data["9"] = [np.mean(data.values[0][239:256]), "NaN", "NaN"]

            mean_stimulated_data = pd.concat([data.iloc[:, i:i + 3].mean(axis=1) for i in range(256, 256 + 9 * 3, 3)],
                                             axis=1)
            mean_stimulated_data["9"] = [np.mean(data.values[0][256:289]), "NaN", "NaN"]

        print("mean baseline")
        print(mean_baseline_data)
        print("mean stim")
        print(mean_stimulated_data)

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

        # if rows containing NaNs were found, delete all rows with this id
        # (if id = 1 in baseline is NaN, corresponding stimulated values are also deleted)
        if rows_with_NaN.empty == False:
            for id in rows_with_NaN["id"]:
                mean_data = mean_data[mean_data["id"] != id]

        # print(mean_data)

        # list of processed data frames
        data_finished.append(mean_data)
    return data_finished, name_list




def export_csv(path_to_folder):
    data_finished, name_list = format_df(path_to_folder)
    print(name_list)
    print(data_finished)
    for i, table in enumerate(data_finished):
        # print(table)
        curve_results_folder = path_to_folder + "Formatted compare data/"
        if not os.path.isdir(curve_results_folder):
            os.makedirs(curve_results_folder)
        table.to_csv(curve_results_folder + name_list[i][0:len(name_list[i]) - 4] + "_format.csv",
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

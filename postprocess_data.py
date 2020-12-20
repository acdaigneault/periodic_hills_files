"""
Name   : postprocess_data.py
Author : Audrey Collard-Daigneault
Date   : 07-12-2020
Desc   : This code plots simulation data with Lethe and other data from literature
         got with Engauge Digitizer or text files.
         (<u>/<ub>, <v>/<ub>, <u'u'>/<ub²>, <v'v'>/<ub²>, <w'w'>/<ub²>, <u'v'>/<ub²>)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

################################################ FILL OUT THIS PART ###################################################

# Reynolds number of the simulation (Currently available for Re = 5600 only)
Re = 5600

# Information about the lethe data
path_to_lethe_data = "/home/audrey/Documents/data_simulation/reynolds_5600/"
file_names_lethe_data = ["data_5_short", "data_3"]  # add all lethe files in this list
# *** NOTE : If writing csv, put file with more rows first because when concaterated columns, all data over the number
#            of rows of the first column won't be added.

# Information about the literature data
path_to_literature_data = "/home/audrey/Documents/data_literature/reynolds_5600/"

# Saving file type ("graph" or "csv")
file_type = "csv"

# Path to save graphs or csv and prefix name
path_to_save = "/home/audrey/Documents/graph/graph/"
name_to_save = "csvfile"

# x/h position with literature data files
x_available = [1, 2, 3, 4, 5, 6, 7, 8]

#######################################################################################################################


def data_to_graph(x_available, Re, path_to_lethe_data, file_names_lethe_data, path_to_literature_data, path_to_save,
                  file_type):
    file_nb = 0

    # Data type to plot
    all_data_type = ["average_velocity_0", "average_velocity_1", "reynolds_normal_stress_0",
                     "reynolds_normal_stress_1", "reynolds_normal_stress_2", "reynolds_shear_stress_uv"]

    # Associate x label to data type
    x_labels = [r"$\langle u \rangle/u_{b}$", r"$\langle v \rangle/u_{b}$", r"$\langle u'u' \rangle/u_{b}^{2}$",
                r"$\langle v'v' \rangle/u_{b}^{2}$", r"$\langle w'w' \rangle/u_{b}^{2}$",
                r"$\langle <u'v' \rangle/u_{b}^{2}$"]

    # Reading Lethe data
    lethe_csv = []
    lethe_data = []
    for lethe_file in file_names_lethe_data:
        lethe_csv.append(path_to_lethe_data + lethe_file + ".csv")
        lethe_data.append(pd.read_csv(lethe_csv[-1], usecols=["Points_0", "Points_1"] + all_data_type, sep=","))

    # A plot for every x value and data type
    for x_value in x_available:
        for index_data_type, data_type in enumerate(all_data_type):
            # Processing Lethe data
            data_to_plot = []
            y_to_plot = []
            for data in lethe_data:
                # Taking data of x value and data type and sorting it for y
                data = data[["Points_0", "Points_1", data_type]]
                y_data = data.loc[(np.abs(data["Points_0"] - x_value) < 1e-6)]

                # If there's no data at x_values, it does fing the 2 nearest x
                nb_unique_x = len(np.unique(y_data["Points_0"]))
                if nb_unique_x is not 1:
                    # Find all x in tolerance = 0.1
                    if nb_unique_x < 1:
                        y_data = data.loc[(np.abs(data["Points_0"] - x_value) < 0.1)]

                    # Get index of value by sorted difference with x_value
                    unique_values = np.unique(y_data["Points_0"])
                    delta = np.abs(unique_values - x_value)
                    index_sorted_delta = np.argsort(delta)

                    # Find the values min and max related to the x_value
                    min_x_value = None
                    max_x_value = None

                    for index in index_sorted_delta:
                        if unique_values[index] < x_value and min_x_value is None:
                            min_x_value = unique_values[index]
                        elif unique_values[index] > x_value and max_x_value is None:
                            max_x_value = unique_values[index]

                    y_data_min = data.loc[(np.abs(data["Points_0"] - min_x_value) < 1e-6)]
                    y_data_max = data.loc[(np.abs(data["Points_0"] - max_x_value) < 1e-6)]
                    y_data = pd.DataFrame(data=y_data_max, columns=y_data_max.columns, index=y_data_max.index)

                    # Linear interpolation : u = (u_2 - u_1) * (x - x_1) / (x2 - x1) + u_1
                    y_data[data_type].sub(y_data_min[data_type])
                    y_data[data_type].mul((x_value - min_x_value) / (max_x_value - min_x_value))
                    y_data[data_type].add(y_data_min[data_type])

                y_data = y_data.sort_values("Points_1")

                # Initialization of arrays prior processing
                data_type_values = np.zeros(1)
                y_values = np.zeros(1)
                nb_z_point = np.zeros(1)  # to count number of points in z axis for each y for the average

                # Initialization of some values prior processing
                last_y = y_data["Points_1"].iloc[0]  # first y point
                y_values[0] = last_y

                # Averaging data in the z axis for each y point
                for index, y in y_data["Points_1"].iteritems():
                    if ~np.isclose(y, last_y):  # next y point
                        data_type_values = np.append(data_type_values, 0)
                        nb_z_point = np.append(nb_z_point, 0)
                        y_values = np.append(y_values, y)

                    data_type_values[-1] += y_data[data_type][index]
                    nb_z_point[-1] += 1
                    last_y = y

                data_type_values /= nb_z_point

                data_to_plot.append(data_type_values)
                y_to_plot.append(y_values)

            # Taking literature data
            [Breuer2009_data, Rapp2009_data] = literature_data_extraction(x_value, data_type, path_to_literature_data,
                                                                          Re)

            file_nb += 1
            if file_type == "graph":
                # Plotting results
                fig, ax = plt.subplots()

                # If there's Lethe data for this x and this data type
                line_type = '-'  # solid line for the first set of data
                for index, name in enumerate(file_names_lethe_data):
                    if data_to_plot[index] is not None:
                        ax.plot(data_to_plot[index], y_to_plot[index], line_type, label='Lethe - ' + name)
                        line_type = '--'  # dashed lines for other Lethe data

                if Breuer2009_data is not None:
                    ax.plot(Breuer2009_data[0], Breuer2009_data[1], '--', color='xkcd:scarlet', label='Breuer2009')

                if Rapp2009_data is not None:
                    ax.plot(Rapp2009_data[0], Rapp2009_data[1], '--', color='xkcd:black', label='Rapp2009')

                ax.set_title(data_type + " at Re = " + str(Re) + " at x = " + str(x_value))
                ax.set_xlabel(x_labels[index_data_type])
                ax.set_ylabel("$y/h$")
                ax.legend()
                fig.savefig(
                    path_to_save + "_" + str(file_nb).rjust(2, '0') + "_" + data_type + "_x_" + str(x_value) + ".png")
                plt.close(fig)
                ax.clear()
            elif file_type == "csv":
                arrays_to_dataframe = pd.DataFrame()
                for i, lethe_file in enumerate(file_names_lethe_data):
                    arrays_to_dataframe.loc[:, data_type + "_" + lethe_file] = pd.Series(data_to_plot[i])
                    arrays_to_dataframe.loc[:, "y/h_" + lethe_file] = pd.Series(y_to_plot[i])

                if Breuer2009_data is not None:
                    arrays_to_dataframe.loc[:, data_type + "_Breuer2009"] = pd.Series(Breuer2009_data[0])
                    arrays_to_dataframe.loc[:, "y/h_Breuer2009"] = pd.Series(Breuer2009_data[1])

                if Rapp2009_data is not None:
                    arrays_to_dataframe.loc[:, data_type + "_Rapp2009"] = pd.Series(Rapp2009_data[0])
                    arrays_to_dataframe.loc[:, "y/h_Rapp2009"] = pd.Series(Rapp2009_data[1])

                arrays_to_dataframe.to_csv(
                    path_to_save + "_" + str(file_nb).rjust(2, '0') + "_" + data_type + "_x_" + str(x_value) + ".csv",
                    index=False, header=True)


# Literature data extraction of files associated with x/h
def literature_data_extraction(x_value, data_type, path_to_literature_data, Re):
    if not Re == 5600:
        assert Re == 5600, "Currently available for Re = 5600 only."

    # Available x value for literature data
    Breuer2009_x_positions = [0.5, 2, 4, 8]
    Rapp2009_x_positions = [0.05, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]

    # Setting file number to x value
    if np.isclose(x_value, 0.05):
        Rapp2009_nb = "01"
    elif np.isclose(x_value, 0.5):
        Breuer2009_nb = ["01", "02", "03", "04", "05", "06"]
        Rapp2009_nb = "02"
    elif x_value == 1:
        Rapp2009_nb = "03"
    elif x_value == 2:
        Breuer2009_nb = ["07", "08", "09", "10", "11", "12"]
        Rapp2009_nb = "04"
    elif x_value == 3:
        Rapp2009_nb = "05"
    elif x_value == 4:
        Breuer2009_nb = ["13", "14", "15", "16", "17", "18"]
        Rapp2009_nb = "06"
    elif x_value == 5:
        Rapp2009_nb = "07"
    elif x_value == 6:
        Rapp2009_nb = "08"
    elif x_value == 7:
        Rapp2009_nb = "09"
    elif x_value == 8:
        Breuer2009_nb = ["19", "20", "21", "22", "23", "24"]
        Rapp2009_nb = "10"

    # Setting file number or column to data type
    if data_type == "average_velocity_0":
        Breuer2009_index = 0
        Rapp2009_column = "u/u_b"
    elif data_type == "average_velocity_1":
        Breuer2009_index = 1
        Rapp2009_column = "v/u_b"
    elif data_type == "reynolds_normal_stress_0":
        Breuer2009_index = 2
        Rapp2009_column = "u'u'/u_b^2"
    elif data_type == "reynolds_normal_stress_1":
        Breuer2009_index = 3
        Rapp2009_column = "v'v'/u_b^2"
    elif data_type == "reynolds_normal_stress_2":
        Breuer2009_index = 4
    elif data_type == "reynolds_shear_stress_uv":
        Breuer2009_index = 5
        Rapp2009_column = "u'v'/u_b^2"

    # Setting False to index variables if not assigned
    if 'Breuer2009_nb' not in locals():
        Breuer2009_nb = None

    if 'Rapp2009_nb' not in locals():
        Rapp2009_nb = None

    if 'Breuer2009_index' not in locals():
        Breuer2009_index = None

    if 'Rapp2009_column' not in locals():
        Rapp2009_column = None

    # Getting Breuer2009 data
    if (x_value in Breuer2009_x_positions) and (Breuer2009_nb is not None) and (Breuer2009_index is not None):
        Breuer2009_nb = Breuer2009_nb[Breuer2009_index]
        Breuer2009_csv = path_to_literature_data + "Breuer2009/Breuer2009_" + str(
            Breuer2009_nb) + ".csv"
        Breuer2009_data = pd.read_csv(Breuer2009_csv, usecols=["x", "Curve" + str(Breuer2009_nb)], sep=",")
        Breuer2009_data = [np.array(Breuer2009_data["Curve" + str(Breuer2009_nb)]), np.array(Breuer2009_data["x"])]
    else:
        Breuer2009_data = None

    # Getting Rapp2009 data
    if x_value in Rapp2009_x_positions and (Rapp2009_column is not None) and (Rapp2009_nb is not None):
        Rapp2009_csv = path_to_literature_data + "Rapp2009/Rapp2009_" + str(Rapp2009_nb) + ".csv"
        Rapp2009_data = pd.read_csv(Rapp2009_csv,
                                    usecols=["y/h", Rapp2009_column], sep=",")
        Rapp2009_data = [np.array(Rapp2009_data[Rapp2009_column]), np.array(Rapp2009_data["y/h"])]
    else:
        Rapp2009_data = None

    return Breuer2009_data, Rapp2009_data


# Data to graph for each x available
data_to_graph(x_available, Re, path_to_lethe_data, file_names_lethe_data, path_to_literature_data,
              path_to_save + name_to_save, file_type)

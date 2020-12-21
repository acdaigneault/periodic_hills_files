"""
Name   : postprocess_data.py
Author : Audrey Collard-Daigneault
Date   : 07-12-2020
Desc   : This code plots simulation data with Lethe and other data from literature
         got with Engauge Digitizer or text files.
         (<u>/<ub>, <v>/<ub>, <u'u'>/<ub²>, <v'v'>/<ub²>, <u'v'>/<ub²>)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

################################################ FILL OUT THIS PART ###################################################

# Reynolds number of the simulation (Currently available for Re = 5600 only)
Re = 5600

# Information about the lethe data
path_to_lethe_data = "C:/Users/Acdai/OneDrive - polymtl.ca/Polytechnique/Session A2020/Periodic Hills Benchmark Case/" \
                     "Data/data_simulation/reynolds_5600/"
file_names_lethe_data = ["data_3"]  # add all lethe files in this list

# Information about the literature data
path_to_literature_data = "C:/Users/Acdai/OneDrive - polymtl.ca/Polytechnique/Session A2020/" \
                          "Periodic Hills Benchmark Case/Data/data_literature/reynolds_5600/"

# Saving file type ("graph" or "csv")
file_type = "csv"

# Path to save graphs or csv and prefix name
path_to_save = "C:/Users/Acdai/OneDrive - polymtl.ca/Polytechnique/Session A2020/Periodic Hills Benchmark Case/Data/" \
               "csv_files_postprocessing/"
name_to_save = "csvwithextra"

# x/h position with literature data files
x_available = [0.05, 1, 2, 3, 4, 5, 6, 7, 8]


#######################################################################################################################


def data_to_graph(x_available, Re, path_to_lethe_data, file_names_lethe_data, path_to_literature_data, path_to_save,
                  file_type):
    file_nb = 0

    # Data type to plot
    all_data_type = ["average_velocity_0", "average_velocity_1", "reynolds_normal_stress_0",
                     "reynolds_normal_stress_1", "reynolds_shear_stress_uv"]


    # Associate x label to data type
    x_labels = [r"$\langle u \rangle/u_{b}$", r"$\langle v \rangle/u_{b}$", r"$\langle u'u' \rangle/u_{b}^{2}$",
                r"$\langle v'v' \rangle/u_{b}^{2}$", r"$\langle <u'v' \rangle/u_{b}^{2}$"]

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
                if nb_unique_x != 1:
                    # Find all x in tolerance = 0.1
                    if nb_unique_x < 1:
                        y_data = data.loc[(np.abs(data["Points_0"] - x_value) < 0.1)]

                    # Get index of value by sorted difference with x_value
                    unique_values = np.unique(y_data["Points_0"])
                    delta = np.abs(unique_values - x_value)
                    index_sorted_delta = np.argsort(delta)

                    # Find the values min and max related to the x_value
                    min_x_value = max_x_value = None

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
                # Find the max value of the array length
                data_size = []
                Breuer2009_data is not None and data_size.append(Breuer2009_data[0].size)
                Rapp2009_data is not None and data_size.append(Rapp2009_data[0].size)
                for i in range(0, len(file_names_lethe_data)):
                    data_size.append(data_to_plot[i].size)

                index_max_size = np.max(np.array(data_size))
                arrays_to_dataframe = pd.DataFrame(index=range(1, index_max_size))

                for i, lethe_file in enumerate(file_names_lethe_data):
                    arrays_to_dataframe.loc[:, data_type + "_" + lethe_file] = pd.Series(data_to_plot[i])
                    arrays_to_dataframe.loc[:, "y_" + lethe_file] = pd.Series(y_to_plot[i])

                if Breuer2009_data is not None:
                    arrays_to_dataframe.loc[:, data_type + "_Breuer2009"] = pd.Series(Breuer2009_data[0])
                    arrays_to_dataframe.loc[:, "y_Breuer2009"] = pd.Series(Breuer2009_data[1])

                if Rapp2009_data is not None:
                    arrays_to_dataframe.loc[:, data_type + "_Rapp2009"] = pd.Series(Rapp2009_data[0])
                    arrays_to_dataframe.loc[:, "y_Rapp2009"] = pd.Series(Rapp2009_data[1])

                arrays_to_dataframe.to_csv(
                    path_to_save + "_" + str(file_nb).rjust(2, '0') + "_" + data_type + "_x_" + str(x_value) + ".csv",
                    index=False, header=True)


# Literature data extraction of files associated with x/h
def literature_data_extraction(x_value, data_type, path_to_literature_data, Re, extra):
    assert Re == 5600, "Currently available for Re = 5600 only."

    # Setting file number to x value
    if np.isclose(x_value, 0.05):
        literature_data_nb = "01"
    elif np.isclose(x_value, 0.5):
        literature_data_nb = "02"
    elif x_value == 1:
        literature_data_nb = "03"
    elif x_value == 2:
        literature_data_nb = "04"
    elif x_value == 3:
        literature_data_nb = "05"
    elif x_value == 4:
        literature_data_nb = "06"
    elif x_value == 5:
        literature_data_nb = "07"
    elif x_value == 6:
        literature_data_nb = "08"
    elif x_value == 7:
        literature_data_nb = "09"
    elif x_value == 8:
        literature_data_nb = "10"
    else:
        literature_data_nb = None

    # Setting file number or column to data type
    if data_type == "average_velocity_0":
        literature_data_type = "u/u_b"
    elif data_type == "average_velocity_1":
        literature_data_type = "v/u_b"
    elif data_type == "reynolds_normal_stress_0":
        literature_data_type = "u'u'/u_b^2"
    elif data_type == "reynolds_normal_stress_1":
        literature_data_type = "v'v'/u_b^2"
    elif data_type == "reynolds_shear_stress_uv":
        literature_data_type = "u'v'/u_b^2"
    else:
        literature_data_type = None


    # Getting literature data
    if literature_data_nb is not None and literature_data_type is not None:
        Rapp2009_csv = path_to_literature_data + "Rapp2009_UFR3-30/Rapp2009_" + str(literature_data_nb) + ".csv"
        Rapp2009_data = pd.read_csv(Rapp2009_csv, usecols=["y/h", literature_data_type], sep=",")
        Rapp2009_data = [np.array(Rapp2009_data[literature_data_type]), np.array(Rapp2009_data["y/h"])]

        Breuer2009_csv = path_to_literature_data + "Breuer2009_UFR3-30/Breuer2009_3-30_" + str(
            literature_data_nb) + ".csv"
        Breuer2009_data = pd.read_csv(Breuer2009_csv, usecols=["y/h", literature_data_type], sep=";")
        Breuer2009_data = [np.array(Breuer2009_data[literature_data_type]), np.array(Breuer2009_data["y/h"])]
    else:
        Rapp2009_data = Breuer2009_data = None

    return Breuer2009_data, Rapp2009_data


# Data to graph for each x available
data_to_graph(x_available, Re, path_to_lethe_data, file_names_lethe_data, path_to_literature_data,
              path_to_save + name_to_save, file_type)

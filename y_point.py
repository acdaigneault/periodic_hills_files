"""
Name   : y_point.py
Author : Audrey Collard-Daigneault
Date   : 16-12-2020
Desc   : This code plots the non-dimensional wall distance (y+) at the lower wall.
         NOTE : yes, this may take a while if your csv files are super duper big :)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

############################# FILL OUT THIS PART #################################
# Reynolds number of the simulation (Currently available for Re = 5600 only)
Re = 5600

# Information about the lethe data
path_to_lethe_data = "/mnt/DATA/niagara_files/"
file_names_lethe_data = ["data_3", "data_5_short"]  # add all lethe files in this list

# Information about the literature data
path_to_literature_data = "/home/audrey/Documents/data_literature/reynolds_5600/"

# Saving file type ("graph" or "csv")
file_type = "graph"

# Path to save graph or csv
path_to_save = "/home/audrey/Documents/graph/graph_3/"
name_to_save = "yplus"

##################################################################################

# Array initiation for files
lethe_csv = []
lethe_data = []
y_plus_per_file = []
x_per_file = []

# Reading Lethe data and sort them by x values
for lethe_file in file_names_lethe_data:
    lethe_csv = path_to_lethe_data + lethe_file + ".csv"
    lethe_data = pd.read_csv(lethe_csv, usecols=["Points_0", "Points_1", "average_velocity_0"], sep=",").sort_values(
        "Points_0")

    u_values = []
    nb_z_point = []
    point_y = []
    point_x = []

    for x in np.unique(lethe_data["Points_0"]):
        x_data = lethe_data.loc[(np.abs(lethe_data["Points_0"] - x) < 1e-6)].sort_values("Points_1")

        count = 0
        y_value = -1
        for i, y in x_data["Points_1"].iteritems():
            if ~np.isclose(y, y_value):
                count += 1
                if count <= 3:
                    u_values.append(0)
                    nb_z_point.append(0)
                    point_x.append(x)
                    point_y.append(y)

            if count <= 3:
                u_values[-1] += x_data["average_velocity_0"][i]
                nb_z_point[-1] += 1

            y_value = y

    u_values = np.array(u_values) / np.array(nb_z_point)
    dudy = []
    viscosity = 1.78571E-04
    y_plus = []
    x = []
    sep_reat = []

    for i in range(0, len(u_values), 3):
        poly = lagrange(np.array(point_y[i:i + 3]), np.array(u_values[i:i + 3]))
        poly_diff = np.polyder(poly)
        dudy.append(poly_diff(point_y[i]))
        dy = point_y[i + 1] - point_y[i]
        y_plus.append(np.sqrt(np.abs(dudy[-1] * viscosity)) * (dy / 2) / viscosity)
        x.append(point_x[i])
        if i != 0 and dudy[-1] * dudy[-2] <= 0:
            sep_reat.append((x[-1] + x[-2]) / 2)

    y_plus_per_file.append(y_plus)
    x_per_file.append(x)

# Literature data from Breuer2009
breuer2009_csv = path_to_literature_data + "Breuer2009/Breuer2009_27.csv"
literature_data = pd.read_csv(breuer2009_csv, usecols=["x", "Curve27"], sep=",")
literature_name = "Breuer2009"

if file_type == "graph":
    fig, ax = plt.subplots()

    for i, lethe_file in enumerate(file_names_lethe_data):
        if lethe_file == "data_3":
            name = "Lethe - ~1M cells - full simulation"
        else:
            name = "Lethe - ~4M cells - half simulation"
        ax.plot(x_per_file[i], y_plus_per_file[i], label=name)
    ax.plot(literature_data["x"], literature_data["Curve27"], '--', color='xkcd:scarlet', label=literature_name)
    ax.set_title("Distribution of $y^{+}$ along the lower wall at Re = " + str(Re))
    ax.set_xlabel("$x/h$")
    ax.set_ylabel("$y^{+}$")
    ax.legend()
    fig.savefig(path_to_save + name_to_save + ".png")
    plt.show()
    plt.close(fig)
    ax.clear()
elif file_type == "csv":
    arrays_to_dataframe = pd.DataFrame(
        {"x_" + literature_name: literature_data["x"], "yplus_" + literature_name: literature_data["Curve27"]})

    for i, lethe_file in enumerate(file_names_lethe_data):
        arrays_to_dataframe.loc[:, 'x_' + lethe_file] = pd.Series(x_per_file[i])
        arrays_to_dataframe.loc[:, 'yplus_' + lethe_file] = pd.Series(y_plus_per_file[i])

    arrays_to_dataframe.to_csv(path_to_save + name_to_save + ".csv", index=False, header=True)

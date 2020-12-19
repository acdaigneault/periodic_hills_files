"""
Name   : postprocess_data.py
Author : Audrey Collard-Daigneault
Date   : 16-12-2020
Desc   : This code pl
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import re

# Periodic Hills lines
# Polynomial coefficients
a1 = 2.800000000000E+01
b1 = 0.000000000000E+00
c1 = 6.775070969851E-03
d1 = -2.124527775800E-03

a2 = 2.507355893131E+01
b2 = 9.754803562315E-01
c2 = -1.016116352781E-01
d2 = 1.889794677828E-03

a3 = 2.579601052357E+01
b3 = 8.206693007457E-01
c3 = -9.055370274339E-02
d3 = 1.626510569859E-03

a4 = 4.046435022819E+01
b4 = -1.379581654948E+00
c4 = 1.945884504128E-02
d4 = -2.070318932190E-04

a5 = 1.792461334664E+01
b5 = 8.743920332081E-01
c5 = -5.567361123058E-02
d5 = 6.277731764683E-04

a6 = 5.639011190988E+01
b6 = -2.010520359035E+00
c6 = 1.644919857549E-02
d6 = 2.674976141766E-05

H = 28.0
max_x = 9 * H
max_y = 3.035 * H
x_vector = np.linspace(0, max_x, 100)
y = np.empty(len(x_vector))

for i, x in enumerate(x_vector):
    new_x = (max_x - x)

    # Polynomial equations :
    if 0 <= x < 9:
        y[i] = a1 + b1 * x + c1 * x ** 2 + d1 * x ** 3
        # Checking if y is under H and correction
        if y[i] > H:
            y[i] = H

    elif 9 <= x < 14:
        y[i] = a2 + b2 * x + c2 * x ** 2 + d2 * x ** 3

    elif 14 <= x < 20:
        y[i] = a3 + b3 * x + c3 * x ** 2 + d3 * x ** 3

    elif 20 <= x < 30:
        y[i] = a4 + b4 * x + c4 * x ** 2 + d4 * x ** 3

    elif 30 <= x < 40:
        y[i] = a5 + b5 * x + c5 * x ** 2 + d5 * x ** 3

    elif 40 <= x < 54:
        y[i] = a6 + b6 * x + c6 * x ** 2 + d6 * x ** 3
        if y[i] < 0:
            y[i] = 0.0

    elif 243 < x <= 252:
        y[i] = a1 + b1 * new_x + c1 * new_x ** 2 + d1 * new_x ** 3
        # Checking if y is under H and correction
        if y[i] > H:
            y[i] = H

    elif 238 < x <= 243:
        y[i] = a2 + b2 * new_x + c2 * new_x ** 2 + d2 * new_x ** 3

    elif 232 < x <= 238:
        y[i] = a3 + b3 * new_x + c3 * new_x ** 2 + d3 * new_x ** 3

    elif 222 < x <= 232:
        y[i] = a4 + b4 * new_x + c4 * new_x ** 2 + d4 * new_x ** 3

    elif 212 < x <= 222:
        y[i] = a5 + b5 * new_x + c5 * new_x ** 2 + d5 * new_x ** 3

    elif 198 < x <= 212:
        y[i] = a6 + b6 * new_x + c6 * new_x ** 2 + d6 * new_x ** 3
        if y[i] < 0:
            y[i] = 0.0

    else:
        y[i] = 0.0

x_vector /= H
y_bottom = y / H
y_top = max_y * np.ones(len(x_vector)) / H

fig, ax = plt.subplots()

ax.plot(x_vector, y_bottom, '-k', linewidth=0.5)
ax.plot(x_vector, y_top, '-k', linewidth=0.5)

# Data type to plot
all_data_type = ["average_velocity_0", "average_velocity_1", "reynolds_normal_stress_0",
                 "reynolds_normal_stress_1", "reynolds_shear_stress_uv"]

# Associate x label to data type
x_labels = [r"$\langle u \rangle/u_{b}$", r"$\langle v \rangle/u_{b}$", r"$\langle u'u' \rangle/u_{b}^{2}$",
            r"$\langle v'v' \rangle/u_{b}^{2}$", r"$\langle u'v' \rangle/u_{b}^{2}$"]

path = "/home/audrey/Documents/graph/graph/"
name_saved = "csvnoBreuer"
range_of_files = list(range(1, 49))
range_of_files = [str(int).rjust(2, '0') for int in range_of_files]

index = 4
name = all_data_type[index]
label = x_labels[index]
factor = 20

filenames = []
bla = False
for file_nb in range_of_files:
    file_name = [filename for filename in os.listdir(path) if filename.startswith(name_saved + "_" + file_nb)][0]
    data_type = file_name.split(name_saved + "_" + file_nb + "_")[1].split("_x")[0]
    x = float(file_name.split("_x_")[1].split(".csv")[0])

    if data_type == name:
        data = pd.read_csv(path + file_name, sep=",")

        for i in range(0, 6, 2):
            data[data.columns[i]] = factor * data[data.columns[i]] + x

        if bla == False:
            ax.plot(data[data.columns[2]], data[data.columns[3]], '--', color='xkcd:crimson',
                    label="Lethe - ~1M cells - full simulation",
                    linewidth=0.75)
            ax.plot(data[data.columns[0]], data[data.columns[1]], '--', color='xkcd:bright blue',
                    label="Lethe - ~4M cells - half simulation",
                    linewidth=0.75)
            ax.plot(data[data.columns[4]], data[data.columns[5]], '--', color='xkcd:black', label="Experimental",
                    linewidth=0.75)
        else:
            ax.plot(data[data.columns[2]], data[data.columns[3]], '--', color='xkcd:crimson',
                    linewidth=0.75)
            ax.plot(data[data.columns[0]], data[data.columns[1]], '--', color='xkcd:bright blue',
                    linewidth=0.75)
            ax.plot(data[data.columns[4]], data[data.columns[5]], '--', color='xkcd:black',
                    linewidth=0.75)
        bla = True

ax.set_title(name + " at Re = 5600")
ax.set_xlabel(str(factor) + "$*$" + label)
ax.set_ylabel("$y/h$")
plt.vlines([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 0, 3.035, linestyle=':', color='xkcd:dark grey', linewidth=0.75)
ax.set_ybound(-0.5, 4.5)
plt.gca().set_aspect('equal', adjustable='box')
ax.legend(fontsize='x-small')
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.savefig(path + "on_geometry/" + name + ".png", dpi=300)
plt.close(fig)
ax.clear()

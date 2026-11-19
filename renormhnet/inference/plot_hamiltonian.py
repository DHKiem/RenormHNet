# plot hamiltonian reading Hamiltonian.csv
# Author: DHKiem


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plotjx():
    lattice = np.loadtxt("lat.dat")
    positions = np.loadtxt("site_positions.dat")
    lattice = lattice.transpose()
    positions = positions.transpose()

    ######### Symmed Jx read  ###########
    sym_jx = pd.read_csv("Hamiltonian.csv")
    ######### Symmed Jx end  ###########


    # calculate distance from Rx, Ry, Rz
    def calculate_distance(row):
        a1 = int(row['Atom1']-1)
        a2 = int(row['Atom2']-1)

        #print("atom1 pos:", atom1position)
        #print("atom2 pos:", atom2position)
        position_vec = row['Rx'] *lattice[0] + row['Ry'] *lattice[1] + row['Rz'] *lattice[2] + (positions[a2] - positions[a1])
        distance = np.linalg.norm(position_vec)
        return distance

    ######### Plotting  ###########

    sym_jx['distance'] = sym_jx.apply(calculate_distance, axis=1)
    sym_jx = sym_jx.sort_values(by='distance')
    sum_jx = sum(sym_jx['J'] )

    plt.plot(sym_jx['distance'], sym_jx['J'], color='royalblue',  marker='o', linestyle='--', alpha=1.0)
    plt.xlabel('Distance (Angstrom)', fontsize=14)
    plt.ylabel('J (meV)', fontsize=14)
    plt.tight_layout()
    plt.savefig("jx_plot.png", dpi=300)
    print("Plot saved as jx_plot.png")
import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import cm
plt.rcParams['text.usetex'] = True
import sys
import time
import scipy.interpolate as interp
import math


def plot_specific_heat(N):
    print("plotting specific heat (3D) ...")

    fig = plt.figure()
    ax = plt.axes(projection='3d')


    for filename in os.listdir("results/3DData/" + N + "/C/"):
        if filename == "dummyFile.txt":
            continue
        J = float(os.path.splitext(filename)[0])
        file = open("results/3DData/" + N + "/C/" + filename, 'r')
        lines = file.readlines()
        X = []
        Y = []
        Z = []
        for line in lines:
            x, z = line.split("\t")
            X += [float(x)]
            Y += [J]
            if math.isnan(float(z)):
                z = 0.0
            Z += [float(z)]
        ax.plot(X, Y, Z, lw = 1, ls = "solid", color = "blue", alpha = 1.0)

    ax.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ für $N$ = ' + N, fontsize = 18)
    ax.set_xlabel(r'$T$ in $J_2$/$k_B$', fontsize = 18)
    ax.set_ylabel(r'$J_1$ / $J_2$', fontsize = 18)
    ax.set_zlabel(r'$C/N$ in $J_2$', fontsize = 18)
    ax.legend(loc = 'best' ,frameon = False, fontsize = 14)

    plt.savefig("results/3DData/" + N + "_specific_heat.png")

def plot_susceptibility(N):
    print("plotting suszeptibility (3D) ...")

    fig = plt.figure()
    ax = plt.axes(projection='3d')


    for filename in os.listdir("results/3DData/" + N + "/X/"):
        if filename == "dummyFile.txt":
            continue
        J = float(os.path.splitext(filename)[0])
        file = open("results/3DData/" + N + "/X/" + filename, 'r')
        lines = file.readlines()
        X = []
        Y = []
        Z = []
        for line in lines:
            x, z = line.split("\t")
            X += [float(x)]
            Y += [J]
            if math.isnan(float(z)):
                z = 0.0
            Z += [float(z)]
        ax.plot(X, Y, Z, lw = 1, ls = "solid", color = "blue", alpha = 1.0)

    ax.set_title(r'suszepibilität $\\chi/N$ für $N$ = ' + N, fontsize = 18)
    ax.set_xlabel(r'$T$ in $J_2$/$k_B$', fontsize = 18)
    ax.set_ylabel(r'$J_1$ / $J_2$', fontsize = 18)
    ax.set_zlabel(r'$\\chi/N$ in $J_2$', fontsize = 18)
    ax.legend(loc = 'best' ,frameon = False, fontsize = 14)

    plt.savefig("results/3DData/" + N + "_specific_heat.png")


if __name__ == "__main__":

    start_time = time.time()

    N = sys.argv[1]

    print("N = " + str(N))
    
    plot_specific_heat(N)
    #plot_susceptibility(N)

    end_time = time.time()

    print("plotting done. This took %.2f seconds" % (end_time - start_time) )

    if (sys.argv[2] == "show"):
        plt.show()
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
        h, J = filename.split("_")
        h = float(h)
        J = float(J.rstrip(".txt"))
        file = open("results/3DData/" + N + "/C/" + filename, 'r')
        lines = file.readlines()
        X = []
        Y = []
        Z = []
        for line in lines:
            x, z = line.split("\t")
            if "nan" in z or "inf" in z:#math.isnan(float(z)):
                continue
            X += [float(x)]
            Y += [J]
            Z += [float(z)]
        ax.plot(X, Y, Z, lw = 1, ls = "solid", color = "blue", alpha = 1.0)

    ax.set_title(r"spezifische Wärmekapazität pro Spin $C/N$" + "\n" + r"für $N$ = " + N + " mit h = " + str(h), fontsize = 18)
    ax.set_xlabel(r'$T$ $k_B$ / $J_2$', fontsize = 18)
    ax.set_ylabel(r'$J_1$ / $J_2$', fontsize = 18)
    ax.set_zlabel(r'$C/N$ in $J_2$', fontsize = 18)
    #ax.legend(loc = 'best' ,frameon = False, fontsize = 14)
    print(h)

    plt.savefig("results/3DData/" + N + "_specific_heat_" + str(h) + ".png")

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
            if "nan" in z or "inf" in z:#math.isnan(float(z)):
                continue
            X += [float(x)]
            Y += [J]
            Z += [float(z)]
        ax.plot(X, Y, Z, lw = 1, ls = "solid", color = "blue", alpha = 1.0)

    ax.set_title('Suszeptibilität $\\chi/N$ für $N$ = ' + N, fontsize = 18)
    ax.set_xlabel(r'$T$ $k_B$ / $J_2$', fontsize = 18)
    ax.set_ylabel(r'$J_1$ / $J_2$', fontsize = 18)
    ax.set_zlabel('$\\chi/N$ in $J_2$', fontsize = 18)
    #ax.legend(loc = 'best' ,frameon = False, fontsize = 14)

    plt.savefig("results/3DData/" + N + "_susceptibility.png")


if __name__ == "__main__":

    start_time = time.time()

    N = sys.argv[1]
  
    if (sys.argv[3] == "C"):
        plot_specific_heat(N)
    if (sys.argv[4] != "noX"):
        plot_susceptibility(N)

    end_time = time.time()

    print("plotting done. This took %.2f seconds" % (end_time - start_time) )

    if (sys.argv[2] == "show"):
        plt.show()
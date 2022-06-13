import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import scipy.optimize
import scipy.stats
plt.rcParams['text.usetex'] = True

N_color = []
N_color_LOW = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta")]#, ("14", "brown"), ("16", "purple")]#, ("18", "tomato")]
N_color_HIGH = [("18", "tomato"), ("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta"), ("28", "brown"), ("30", "purple"), ("32", "tomato")]

peak_offset = 500 # 1000

# N_color = [("12", "orange")]

def sort_data(X, Y):
    length = len(X)
    for i in range(length-1):
        for j in range(0, length-i-1):
            x1 = X[j]; x2 = X[j+1]
            if x1 > x2:
                X[j], X[j+1] = X[j+1], X[j]
                Y[j], Y[j+1] = Y[j+1], Y[j]
    return X, Y

def expFunc(x: float, A: float, k: float, x_0: float, y_0: float) -> float:
    return A * np.exp(k * (x + x_0)) + y_0

def get_spin_gap(N, J, filename) -> (float, float, float, float):
    file = open("results/" + N + "/data/spin_gap_data/" + filename, 'r')
    lines = file.readlines()
    fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
    X = []; Y = []
    for i in range(5,len(lines)):
        arr = lines[i].split("\t")
        x = arr[0]
        y = arr[1]
        X += [float(x)] # beta -> T
        Y += [float(y)]

    # print(len(X))
    # print(len(Y))

    X.reverse(); Y.reverse()
    X_fit = []; Y_fit = []
    for i in range(0, len(X)-peak_offset):
        # print("%i, %f, %f, %f" %(i, X[i], Y[i+peak_offset], Y[i]))
        # if X[i] < 10**(-5): continue
        if Y[i] > Y[i+peak_offset]:
            #print("break")
            break
        X_fit += [X[i]]; Y_fit += [np.log(Y[i])] # log = ln [np.log(Y[i])]
    if (len(X_fit) == 0):
        for i in range(0, len(X) - int(len(X)/3)):
            X_fit += [X[i]]; Y_fit += [np.log(Y[i])]
    # print(len(X_fit))
    # print(len(Y_fit))
    #X_fit.reverse(); Y_fit.reverse()
    # defaultParams = (1.0, 0.5, -0.05, 0.01)
    # params, cv = scipy.optimize.curve_fit(expFunc, X_fit, Y_fit, defaultParams)
    # A, k, x_0, y_0 = params
    # res = scipy.stats.linregress(X_fit, Y_fit)
    # m = res.slope; b = res.intercept
    m, b = np.polyfit(X_fit, Y_fit, 1)
    subfig2.plot(X, Y, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "blue", label = "QT data")
    # Y_fitted = [expFunc(x, A, k, x_0, y_0) for x in X_fit]
    Y_fitted = [np.exp(m * x + b) for x in X_fit]
    subfig2.plot(X_fit, Y_fitted, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "red", label = "exp fit")
    #subfig2.set_xlabel(r'T in $k_B$ / $J_2$', fontsize = 25)
    subfig2.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
    subfig2.set_ylabel('$\\chi/N$ in $J_2$', fontsize = 25)
    subfig2.set_title("$\\chi/N$ für N = " + N + r" mit $J_1$ / $J_2$ = " + J, fontsize = 25)

    # subfig2.axhline(0, color = "grey")
    subfig2.legend(loc = 'best' ,frameon = False, fontsize = 20)

    plt.axvline(x=X_fit[len(X_fit)-1], color='black', linestyle='--')
    #subfig2.set_xscale("log")
    subfig2.set_yscale("log")
    # plt.xscale("log")
    # plt.yscale('log')
    # plt.xlim(X[0] - 0.005, X_fit[len(X_fit)-1] + 0.005)
    plt.savefig("results/" + N + "/data/spin_gap_data/X_J" + J + ".png")
    plt.close(fig2)
    return np.exp(b), abs(m), 0.0, 0.0
    # return A, k, x_0, y_0
    # return 0.0, 0.0, 0.0, 0.0


if __name__ == "__main__":

    regime = sys.argv[1]
    if regime == "low": N_color = N_color_LOW
    elif regime == "high": N_color = N_color_HIGH
    else: exit()

    print("plotting spin gap ...")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

    used_N = "N"

    for N, c in N_color:
        print("N = " + N)
        # QT results
        lbl = "QT: N = " + N
        X = []; Y = []
        for filename in os.listdir("results/" + N + "/data/spin_gap_data/"):
            if filename == "dummyFile.txt" or filename == "data_placeholder.txt": continue
            if filename[len(filename)-4:] == ".png": continue
            J = filename[len("X_J"):-len(".txt")]
            # print(J)
            A, k, x_0, y_0 = get_spin_gap(N, J, filename)
            # print(str(A) + " " + str(k)  + " " + str(x_0)  + " " + str(y_0))
            X += [float(J)]
            Y += [float(k)]
        X, Y = sort_data(X, Y)
        subfig1.plot(X, Y, lw = 1, ls = "dashed", markersize = 2, marker = "o", color = c, label = lbl)
        # ED results
        file = open("results/" + N + "/data/data_spin_gap.txt", 'r') # _data_spin_gap / _data_spin_gap_with_index
        lines = file.readlines()
        lbl = "ED: N = " + N
        X = []; Y = []
        for i in range(7,len(lines)):
            arr = lines[i].split("\t")
            x = arr[0]
            y = arr[1]
            X += [float(x)]
            Y += [float(y)]
        subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = c, label = lbl, alpha = 0.5)

        used_N += "_" + N

        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 25)
        subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$', fontsize = 25)
        # subfig1.set_title(r'Spingap Energies $\Delta E_{gap}$ für $\Delta = 1$', fontsize = 18)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

        plt.savefig("results/" + "spin_gap_with_QT_" + used_N + ".png")
        #plt.show()

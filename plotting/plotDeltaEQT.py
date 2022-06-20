import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
import os
import sys
import numpy as np
import scipy.optimize
from typing import Tuple
import gc
plt.rcParams['text.usetex'] = True

N_color = []
N_color_LOW = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color_HIGH = [("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta")]#, ("28", "brown"), ("30", "purple"), ("32", "tomato")]

peak_offset = 2000 #1500 # 1000
fit_samples = 5 # 5
search_start_percent = 4/5
search_end_percent = 1/5
max_n = 5 # min = 1; max = 5
no_ED = False

SAVE_FULL_PLOTS = False

def sort_data(X, Y, A):
    length = len(X)
    for i in range(length-1):
        for j in range(0, length-i-1):
            x1 = X[j]; x2 = X[j+1]
            if x1 > x2:
                X[j], X[j+1] = X[j+1], X[j]
                Y[j], Y[j+1] = Y[j+1], Y[j]
                A[j], A[j+1] = A[j+1], A[j]
    return X, Y, A

def expFunc(x: float, A: float, k: float) -> float:
    return x ** 2 * A * np.exp(k * x)

def get_excitation_energie(n, N, J, filename) -> Tuple[float, float]:
    file = open("results/" + N + "/data/excitation_energies_data/" + str(n+1) + "/" + filename, 'r')
    ED_QT = filename[len(filename)-6:-4]
    lines = file.readlines()
    linesh = "0.0"
    X = []; Y = []
    for i in range(5,len(lines)):
        arr = lines[i].split("\t")
        x = arr[0]
        y = arr[1]
        X += [float(x)]
        Y += [float(y)]
    file.close()
    X.reverse(); Y.reverse()

    X_fit = X[0:int(len(X)*search_start_percent)]; Y_fit = Y[0:int(len(X)*search_start_percent)]
    X_fit = np.array(X_fit); Y_fit = np.array(Y_fit)
    params, cv = scipy.optimize.curve_fit(expFunc, X_fit, Y_fit, (1.0, 0.1))
    A, k = params
    # calc R2
    residuals = Y_fit - expFunc(X_fit, A, k)
    SSE = np.sum(residuals**2)
    diff = Y_fit - Y_fit.mean()
    square_diff = diff ** 2
    SST = square_diff.sum()
    R2 = 1 - SSE/SST
    # set ranges
    X_fit_range = np.array(X_fit); Y_fit_range = np.array(Y_fit)

    for p in range(1, fit_samples+1):
        percent = search_start_percent + (search_end_percent - search_start_percent) * float(p) / float(fit_samples)
        #print(percent)
        X_fit = []; Y_fit = []
        for i in range(0, int(len(X)*percent)):
            X_fit += [X[i]]; Y_fit += [Y[i]] 
        X_fit = np.array(X_fit); Y_fit = np.array(Y_fit)
        params, cv = scipy.optimize.curve_fit(expFunc, X_fit, Y_fit, (1.0, 0.1))
        A_new, k_new = params
        # calc R2
        residuals = Y_fit - expFunc(X_fit, A_new, k_new)
        SSE = np.sum(residuals**2)
        diff = Y_fit - Y_fit.mean()
        square_diff = diff ** 2
        SST = square_diff.sum()
        R2_new = 1 - SSE/SST
        # compare to current best
        if R2_new >= R2:
            R2 = R2_new
            A = A_new
            k = k_new
            X_fit_range = X_fit; Y_fit_range = Y_fit
            Y_fitted = [expFunc(x, A, k) for x in X_fit_range]

    if SAVE_FULL_PLOTS:
        fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))

        subfig2.plot(X_fit_range, Y_fit_range, lw = 1, ls = "solid", markersize = 5, marker = "o", color = "green", label = "range")

        subfig2.plot(X, Y, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "blue", label = "QT data")
        Y_fitted = [expFunc(x, A, k) for x in X_fit_range]
        subfig2.plot(X_fit_range, Y_fitted, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "red", label = "exp fit, R = " + str(R2))
        subfig2.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig2.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 25)
        subfig2.set_title(r'Anregungsenergieren $\Delta E$' + linesh, fontsize = 25)

        # subfig2.axhline(0, color = "grey")
        subfig2.legend(loc = 'best' ,frameon = False, fontsize = 20)

        plt.axvline(x=X_fit_range[len(X_fit_range)-1], color='black', linestyle='--')
        #subfig2.set_xscale("log")
        subfig2.set_yscale("log")
        fig2.savefig("results/" + N + "/data/excitation_energies_data/" + str(n+1) + "/" + "C_J" + J + "_" + ED_QT + ".png")
        plt.cla()
        plt.clf()
        plt.close(fig2)
        del fig2, subfig2
        gc.collect()

    return np.exp(A), abs(k)

def save_excitation_energy_data(N, X, Y, YErr, A, AErr) -> None:
    cnt = -1
    for filename in os.listdir("results/" + N + "/data/"):
        if "excitation_energy_data_" in filename and "AMP" not in filename and "Amp" not in filename:
            _cnt = filename[len("excitation_energy_data_"):-len("_QT.txt")]
            if int(_cnt) > cnt: cnt = int(_cnt)
    outFile = open("results/" + N + "/data/excitation_energy_data_" + str(cnt+1) + "_QT.txt", "w")
    for i in range(len(X)):
        outFile.write("%f\t%ft%f" % (X[i], Y[i], YErr[i]) )
    outFile.close()
    outFile = open("results/" + N + "/data/excitation_energy_data_AMP_" + str(cnt+1) + "_QT.txt", "w")
    for i in range(len(X)):
        outFile.write("%f\t%ft%f" % (X[i], A[i], AErr[i]) )
    outFile.close()
    # save excitation energies
    if SAVE_FULL_PLOTS:
        fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
        subfig3.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = "red", label = lbl)
        X = np.asarray(X)
        Y = np.asarray(Y)
        YErr = np.asarray(YErr)
        subfig3.fill_between(X, Y - YErr, Y + YErr, color = "blue", alpha = 0.2)
        subfig3.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig3.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 25)
        subfig3.set_title(r'Anregungsenergieren $\Delta E$', fontsize = 25)
        subfig3.axhline(0, color = "grey")
        subfig3.legend(loc = 'best' ,frameon = False, fontsize = 20)
        fig3.savefig("results/" + N + "/excitation_energy_data_" + str(cnt+1) + "_QT.png")
        plt.cla()
        plt.clf()
        plt.close(fig3)
        del fig3, subfig3
        # save amplitude
        fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
        subfig3.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = "red", label = lbl)
        X = np.asarray(X)
        Y = np.asarray(Y)
        YErr = np.asarray(YErr)
        subfig3.fill_between(X, Y - YErr, Y + YErr, color = "blue", alpha = 0.2)
        subfig3.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig3.set_ylabel(r'$A$  in $J_2$', fontsize = 25)
        subfig3.set_title(r'Amplituden $A$', fontsize = 25)
        subfig3.axhline(0, color = "grey")
        subfig3.legend(loc = 'best' ,frameon = False, fontsize = 20)
        fig3.savefig("results/" + N + "/excitation_energy_data_AMP_" + str(cnt+1) + "_QT.png")
        plt.cla()
        plt.clf()
        plt.close(fig3)
        del fig3, subfig3
        gc.collect()


if __name__ == "__main__":

    if len(sys.argv) > 1:
        regime = sys.argv[1]
        if regime == "low": N_color = N_color_LOW#; print("low regime")
        elif regime == "high": N_color = N_color_HIGH; no_ED = True#; print("high regime")
        else: N_color = N_color_LOW; print("default low (wrong args)")
    else: N_color = N_color_LOW; print("default low (no args)")

    for i in range(len(N_color)):
        print("plotting excitation energies ...")
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        figAmp, subfigAmp = plt.subplots(1,1,figsize=(16,9))
        used_N = "N"
        for j in range(i, len(N_color)):
            N, c = N_color [j]
    #for N, c in N_color:
            print("N = " + N + ":") 
            # QT results
            print("QT (exp fit)...")
            X_arr = [[]] * max_n; Y_arr = [[]] * max_n; A_arr = [[]] * max_n
            for n in range(0, max_n):
                print("n = " + str(n+1))
                #lbl = "QT (exp fit): N = " + N
                lbl = "N = " + N
                X_arr[n] = []; Y_arr[n] = [];  A_arr[n] = []
                for filename in os.listdir("results/" + N + "/data/excitation_energies_data/" + str(n+1) + "/"):
                    if filename[len(filename)-6:] != "QT.txt": continue
                    J = filename[len("C_J"):-len("QT.txt")]
                    linesh = "0.0"
                    # print(J)
                    A, k = get_excitation_energie(n, N, J, filename)
                    # print(str(A) + " " + str(k)  + " " + str(x_0)  + " " + str(y_0))
                    X_arr[n] += [float(J)]
                    Y_arr[n] += [float(k)]
                    A_arr[n] += [float(A)]
                X_arr[n], Y_arr[n], A_arr[n] = sort_data(X_arr[n], Y_arr[n], A_arr[n])

            X = []; Y = []; YErr = []; A = []; AErr = []
            for pos in range(len(X_arr[0])):
                y = []; a = []
                for n in range(max_n):
                    y += [Y_arr[n][pos]]
                    a += [A_arr[n][pos]]
                y = np.array(y); a = np.array(a)
                X += [X_arr[0][pos]]; Y += [y.mean()]; YErr += [y.std()]; A += [a.mean()]; AErr += [a.std()/a.mean()]

            save_excitation_energy_data(N, X, Y, YErr, A, AErr)

            subfig1.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = c, label = lbl)
            X = np.asarray(X)
            Y = np.asarray(Y)
            YErr = np.asarray(YErr)
            subfig1.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.2)

            subfigAmp.plot(X, AErr, lw = 1, ls = "solid", markersize = 0, marker = "o", color = c, label = lbl, alpha = 0.5)
            # A = np.asarray(A)
            # AErr = np.asarray(AErr)
            # subfigAmp.fill_between(X, A - AErr, A + AErr, color = c, alpha = 0.2)

            if not no_ED:
                # ED results exp fit
                print("ED (exp fit)...")
                lbl = "ED (exp fit): N = " + N
                X = []; Y = []; A_arr = []
                for filename in os.listdir("results/" + N + "/data/excitation_energies_data/1/"):
                    if filename[len(filename)-6:] != "ED.txt": continue
                    J = filename[len("X_J"):-len("ED.txt")]
                    # print(J)
                    A, k = get_excitation_energie(0, N, J, filename)
                    # print(str(A) + " " + str(k)  + " " + str(x_0)  + " " + str(y_0))
                    X += [float(J)]
                    A_arr += [float(A)]
                X, Y, A_arr = sort_data(X, Y, A_arr)
                subfig1.plot(X, Y, lw = 0, ls = "dotted", markersize = 2, marker = "o", color = c, alpha = 0.5)#, label = lbl, alpha = 0.5)
                # subfigAmp.plot(X, A_arr, lw = 0, ls = "dotted", markersize = 2, marker = "o", color = c, alpha = 0.5)
                # ED results
                print("ED (dispersion)...")
                lbl = "ED : N = " + N
                X = []; Y = []
                file = open("results/" + N + "/data/data_delta_E.txt", 'r')
                lines = file.readlines()
                lbl = "ED: N = " + N
                X = []; Y = []
                for i in range(8,len(lines)):
                    arr = lines[i].split("\t")
                    x = arr[0]
                    y = arr[1]
                    X += [float(x)]
                    Y += [float(y)]
                file.close()
                subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 0, marker = "o", color = c, alpha = 0.4) #  label = lbl,
            
            print()

            used_N += "_" + N

            subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
            subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 25)
            subfig1.set_title(r'Anregungsenergieren $\Delta E$ mit '  + str(max_n) + " Startvektoren (QT), h = " + linesh, fontsize = 25)

            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

            fig1.savefig("results/" + "excitation_energies_data_with_QT_" + used_N + ".png")

            subfigAmp.set_yscale("log")

            subfigAmp.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
            subfigAmp.set_ylabel(r'$\sigma_{rel}(A)$ in $J_2$', fontsize = 25)
            subfigAmp.set_title(r'rel. Standardabweichung Amplituden $\sigma(A)$ der exp. Fits mit ' + str(max_n) + " Startvektoren (QT)", fontsize = 25)

            subfigAmp.axhline(0, color = "grey")
            subfigAmp.legend(loc = 'best' ,frameon = False, fontsize = 20)

            figAmp.savefig("results/" + "spin_gap_with_QT_AMP_" + used_N + ".png")

            # gc.collect(generation=2)

            #plt.show()

        plt.close(fig1)
        plt.close(figAmp)
        # plt.close('all')
        # plt.cla()
        # plt.clf()
        # del fig1, subfig1
        # del figAmp, subfigAmp
        # gc.collect()
        exit()

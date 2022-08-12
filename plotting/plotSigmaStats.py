import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import sys
import os
import gc
from typing import Tuple
import scipy.optimize
plt.rcParams['text.usetex'] = True

import matplotlib as mpl
mpl.rcParams['figure.max_open_warning'] = 0

N_color = []
N_color_LOW = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color_HIGH = [("18", "tomato"), ("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta"), ("28", "brown"), ("30", "purple"), ("32", "tomato")]
n_color = [("1", "red"), ("2", "blue"), ("3", "green"), ("4", "tomato")]
colors = ["red", "blue", "green", "magenta", "tomato", "brown", "purple"]#, "cyan"]
N_Array = []

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

line_width = 4
marker_size = 0
alph = 0.1

peak_offset = 2000
fit_samples = 3 # 5
search_start_percent = 4/5
search_end_percent = 1/5
max_n = 30

ABS = True
REL = True

MAKE_NEW = False
J_SUM_SCALE = True

valid_step_sizes = [0.1, 0.01, 0.001]
start = 0.0
ends = [0.1, 0.15, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0]
valid_J = [0.010000, 0.358600, 1.006000, 1.753000]

def expFunc(x: float, A: float, k: float) -> float:
    return x * A * np.exp(k * x)

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

def get_spin_gap(n: int, N: int, J: str, filename: str) -> Tuple[float, float]:
    file = open("results/" + N + "/data/spin_gap_data/" + str(n) + "/" + filename, 'r')
    ED_QT = filename[len(filename)-6:-4]
    lines = file.readlines()
    X = []; Y = []
    for i in range(5,len(lines)):
        arr = lines[i].split("\t")
        x = arr[0]
        y = arr[1]
        X += [float(x)] # beta -> T
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
    X_fit_range = np.array(X_fit); Y_fit_range = np.array(Y_fit)

    for p in range(1, fit_samples+1):
        percent = search_start_percent + (search_end_percent - search_start_percent) * float(p) / float(fit_samples)
        X_fit = []; Y_fit = []
        for i in range(0, int(len(X)*percent)):
            X_fit += [X[i]]; Y_fit += [Y[i]] # log = ln [np.log(Y[i])]
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
        if R2_new > R2:
            R2 = R2_new
            A = A_new
            k = k_new
            X_fit_range = X_fit; Y_fit_range = Y_fit

        gc.collect()

    return float(A), abs(k)

def getDataQT(N: str, n: int):
    if MAKE_NEW:
        X = []; Y = [];  A = []
        X, Y, A = newQTdata(N, n)
    else:
        try:
            file = open("results/" + N + "/data/spin_gap_data/" + str(n) + ".txt", 'r')
            lines = file.readlines()
            X = []; Y = [];  A = []
            for i in range(len(lines)):
                x, y, a = lines[i].split("\t")
                X += [float(x)]
                Y += [float(y)]
                A += [float(a)]
            file.close()
        except:
            X = []; Y = [];  A = []
            X, Y, A = newQTdata(N, n)
    return X, Y, A

def newQTdata(N: str, n: int):
    X = []; Y = [];  A = []
    for filename in os.listdir("results/" + N + "/data/spin_gap_data/" + str(n) + "/"):
        if filename[len(filename)-6:] != "QT.txt": continue
        J = filename[len("X_J"):-len("QT.txt")]
        amp, k = get_spin_gap(n, N, J, filename)
        X += [float(J)]
        Y += [float(k)]
        A += [float(amp)]
    X, Y, A = sort_data(X, Y, A)
    file = open("results/" + N + "/data/spin_gap_data/" + str(n) + ".txt", "w")
    for i in range(len(X)):
        file.write(f"{X[i]}\t{Y[i]}\t{A[i]}\n")
    file.close()

def getRel(XDat, YDat, YErrDat):
    X = []; Y = []; YErr = []
    for i in range(len(XDat)):
        if YDat[i] >= 0.000001:
            X += [XDat[i]]; Y += [YDat[i]]; YErr += [YErrDat[i]]
    X = np.asarray(X); Y = np.asarray(Y); YErr = np.asarray(YErr)
    YErr = YErr / Y
    return X, Y, YErr

def avgData(samp, inputdata):
    data_to_avg = []
    for i in range(0, max_n - samp + 1, samp):
        # print("%i, %i, %i, %i" % (max_n, len(inputdata), samp, i))
        data_raw = []
        for j in range(i, i + samp):
            data_raw += [inputdata[j]]
        Y_data = []
        for k in range(0, len(data_raw[0])):
            y = []
            for Yd in data_raw:
                y += [Yd[k]]
            y = np.asarray(y)
            Y_data += [y.mean()]
        data_to_avg += [Y_data]
    Y = []; YErr = []
    for k in range(0, len(data_to_avg[0])):
        y = []
        for Yd in data_to_avg:
            y += [Yd[k]]
        y = np.asarray(y)
        # print("%f, %f" % (y.mean(), y.std()))
        Y += [y.mean()]; YErr +=[y.std()]
    return Y, YErr

def isValid(J):
    for j in valid_J:
        if abs(J-j) <= 0.01:
            return True
    return False

def Npos(N):
    for i in range(len(N_Array)):
        if N_Array[i] == N:
            return i
    return -1

def C_plot_n_for_each_N_sigma(): 
    print("C: plotting dependance of n for fixed N ...")
    abs_rel = ""
    if ABS: abs_rel += "abs"
    if REL: abs_rel += "rel"
    s = -1
    # try:
    for N in N_Array:
        N = str(N)
        # QT results
        J_cnt = 0
        for filename in os.listdir("results/" + N + "/data/excitation_energies_data/1/"):
            x_min = 0
            if ".png" in filename or ".pdf" in filename: continue
            if "ED" in filename: continue
            if "placeholder" in filename: continue
            J = filename[len("C_J"):-len("QT.txt")]
            if not isValid(float(J)): continue
            J_cnt += 1
            print("get data: N: %s, J: %s, %i/%i          \r" % (N, str(J), J_cnt, len(valid_J)), end = "")
            X_arr = [[]] * max_n; Y_arr = [[]] * max_n
            for n in range(max_n):
                file = open("results/" + N + "/data/excitation_energies_data/" + str(n+1) + "/" + filename, 'r')
                lines = file.readlines()
                X = []; X_arr[n] = []
                Y = []; Y_arr[n] = []
                YErr = []
                for i in range(7,len(lines)):
                    x, y = lines[i].split("\t")
                    if float(x) <= 0: continue
                    # if not float(y) > 0.0: continue
                    X += [1.0/float(x)]; X_arr[n] += [1.0/float(x)]
                    Y += [float(y)]; Y_arr[n] += [float(y)]
                file.close()
                if min(X) > x_min: x_min = min(X)
            
            X = X_arr[0]
            filename = "N_" + N + "_n"
            color_count = 0
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            legend = []
            for samp in range(1, int(max_n / 2) + 1):
                s = samp
                print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                Y, YErr = avgData(samp, Y_arr)
                Y = np.asarray(Y); YErr = np.asarray(YErr)
                # ploting
                if ABS:
                    subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                if REL: 
                    Xrel, Yrel, YErrrel = getRel(X, Y, YErr)
                    if min(Xrel) > x_min: x_min = min(Xrel)
                    subfig1.plot(Xrel, YErrrel, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
                legend += [Line2D([0], [0], label = "n = %s" % str(samp), color = colors[color_count], ls = "solid", lw = line_width)]
                color_count += 1
                # saving
                handles, labels = plt.gca().get_legend_handles_labels()
                handles.extend(legend)
                subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
                subfig1.set_ylabel(r'$\sigma$ / $J_2$', fontsize = labelfontsize)
                subfig1.set_title(r"$\sigma$ von $C/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
                subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
                filename += "_" + str(samp)
                if color_count >= len(colors): break
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/sigma/C_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/sigma/C_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...                              ")

def C_plot_N_for_each_n_sigma():
    print("C: plotting dependance of N for fixed n ...")
    abs_rel = ""
    if ABS: abs_rel += "abs"
    if REL: abs_rel += "rel"
    s = -1
    # try:
    for filename in os.listdir("results/" + str(N_Array[0]) + "/data/excitation_energies_data/1/"):
        x_min = 0
        if ".png" in filename or ".pdf" in filename: continue
        if "ED" in filename: continue
        if "placeholder" in filename: continue
        J = filename[len("C_J"):-len("QT.txt")]
        if not isValid(float(J)): continue
        N_SAMP = []
        for N in N_Array:
            # QT results
            J_cnt = 0
            for filename in os.listdir("results/" + str(N) + "/data/excitation_energies_data/1/"):
                if ".png" in filename or ".pdf" in filename: continue
                if "ED" in filename: continue
                if "placeholder" in filename: continue
                Jnew = filename[len("C_J"):-len("QT.txt")]
                if str(J) != str(Jnew): continue
                J_cnt += 1
                # print("get data: N: %s, J: %s, %i/%i          \r" % (N, str(J), J_cnt, len(valid_J)), end = "")
                X_arr = [[]] * max_n; Y_arr = [[]] * max_n
                for n in range(max_n):
                    file = open("results/" + str(N) + "/data/excitation_energies_data/" + str(n+1) + "/" + filename, 'r')
                    lines = file.readlines()
                    X = []; X_arr[n] = []
                    Y = []; Y_arr[n] = []
                    YErr = []
                    for i in range(7,len(lines)):
                        x, y = lines[i].split("\t")
                        if float(x) <= 0: continue
                        # if not float(y) > 0.0: continue
                        X += [1.0/float(x)]; X_arr[n] += [1.0/float(x)]
                        Y += [float(y)]; Y_arr[n] += [float(y)]
                    file.close()
                    if min(X) > x_min: x_min = min(X)
                
                X = X_arr[0]
                color_count = 0
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                sampArr = []
                for samp in range(1, int(max_n / 2) + 1):
                    print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                    Y, YErr = avgData(samp, Y_arr)
                    sampArr += [YErr]
                N_SAMP += [sampArr]
        for samp in range(len(N_SAMP[0])):
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,10))
            legend = []
            filename = "n_" + str(samp+1) + "_N"
            color_count = 0
            for N in range(len(N_Array)):
                print("plotting: %i sample: N: %i, J: %s, %i/%i          \r" % (samp+1, N_Array[N], str(J), J_cnt, len(valid_J)), end = "")
                YErr = N_SAMP[N][samp]
                # ploting
                if ABS:
                    subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                if REL: 
                    Xrel, Yrel, YErrrel = getRel(X, Y, YErr)
                    if min(Xrel) > x_min: x_min = min(Xrel)
                    subfig1.plot(Xrel, YErrrel, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
                legend += [Line2D([0], [0], label = "N = " + str(N_Array[N]), color = colors[color_count], ls = "solid", lw = line_width)]
                color_count += 1
                filename += "_" + str(N_Array[N])
                if color_count >= len(colors): break
            # saving
            handles, labels = plt.gca().get_legend_handles_labels()
            handles.extend(legend)
            subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_ylabel(r'$\sigma$ / $J_2$', fontsize = labelfontsize)
            vec = "einem Startvektor"
            if int(samp) > 0: vec = str(samp+1) + " Startvektoren"
            subfig1.set_title(r"$\sigma$ von $C/N$ bei der QT mit Mittelung" +  "\nüber " + vec + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
            subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/sigma/C_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/sigma/C_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...                                        ")

def X_plot_n_for_each_N_sigma(): 
    print("X: plotting dependance of n for fixed N ...")
    abs_rel = ""
    if ABS: abs_rel += "abs"
    if REL: abs_rel += "rel"
    s = -1
    # try: 
    for N in N_Array:
        N = str(N)
        # QT results
        J_cnt = 0
        for filename in os.listdir("results/" + N + "/data/spin_gap_data/1/"):
            x_min = 0
            if ".png" in filename or ".pdf" in filename: continue
            if "ED" in filename: continue
            if "placeholder" in filename: continue
            J = filename[len("C_J"):-len("QT.txt")]
            if not isValid(float(J)): continue
            J_cnt += 1
            print("get data: N: %s, J: %s, %i/%i          \r" % (N, str(J), J_cnt, len(valid_J)), end = "")
            X_arr = [[]] * max_n; Y_arr = [[]] * max_n
            for n in range(max_n):
                file = open("results/" + N + "/data/spin_gap_data/" + str(n+1) + "/" + filename, 'r')
                lines = file.readlines()
                X = []; X_arr[n] = []
                Y = []; Y_arr[n] = []
                YErr = []
                for i in range(7,len(lines)):
                    x, y = lines[i].split("\t")
                    if float(x) <= 0: continue
                    # if not float(y) > 0.0: continue
                    X += [1.0/float(x)]; X_arr[n] += [1/float(x)]
                    Y += [float(y)]; Y_arr[n] += [float(y)]
                file.close()
                if min(X) > x_min: x_min = min(X)
            
            X = X_arr[0]
            filename = "N_" + N + "_n"
            color_count = 0
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            legend = []
            for samp in range(1, int(max_n / 2) + 1):
                s = samp
                print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                Y, YErr = avgData(samp, Y_arr)
                Y = np.asarray(Y); YErr = np.asarray(YErr)
                # ploting
                if ABS:
                    subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                if REL: 
                    Xrel, Yrel, YErrrel = getRel(X, Y, YErr)
                    if min(Xrel) > x_min: x_min = min(Xrel)
                    subfig1.plot(Xrel, YErrrel, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
                legend += [Line2D([0], [0], label = "n = %s" % str(samp), color = colors[color_count], ls = "solid", lw = line_width)]
                color_count += 1
                # saving
                handles, labels = plt.gca().get_legend_handles_labels()
                handles.extend(legend)
                subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
                subfig1.set_ylabel(r'$\sigma$ / $J_2$', fontsize = labelfontsize)
                subfig1.set_title(r"$\sigma$" + " von $\\chi/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
                subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
                filename += "_" + str(samp)
                if color_count >= len(colors): break
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/sigma/X_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/sigma/X_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...                              ")

def X_plot_N_for_each_n_sigma():
    print("X: plotting dependance of N for fixed n ...")
    abs_rel = ""
    if ABS: abs_rel += "abs"
    if REL: abs_rel += "rel"
    s = -1
    # try:
    for filename in os.listdir("results/" + str(N_Array[0]) + "/data/spin_gap_data/1/"):
        x_min = 0
        if ".png" in filename or ".pdf" in filename: continue
        if "ED" in filename: continue
        if "placeholder" in filename: continue
        J = filename[len("C_J"):-len("QT.txt")]
        if not isValid(float(J)): continue
        N_SAMP = []
        for N in N_Array:
            # QT results
            J_cnt = 0
            for filename in os.listdir("results/" + str(N) + "/data/spin_gap_data/1/"):
                x_min = 0
                if ".png" in filename or ".pdf" in filename: continue
                if "ED" in filename: continue
                if "placeholder" in filename: continue
                Jnew = filename[len("C_J"):-len("QT.txt")]
                if str(J) != str(Jnew): continue
                J_cnt += 1
                print("get data: N: %s, J: %s, %i/%i          \r" % (N, str(J), J_cnt, len(valid_J)), end = "")
                X_arr = [[]] * max_n; Y_arr = [[]] * max_n
                for n in range(max_n):
                    file = open("results/" + str(N) + "/data/spin_gap_data/" + str(n+1) + "/" + filename, 'r')
                    lines = file.readlines()
                    X = []; X_arr[n] = []
                    Y = []; Y_arr[n] = []
                    YErr = []
                    for i in range(7,len(lines)):
                        x, y = lines[i].split("\t")
                        if float(x) <= 0: continue
                        # if not float(y) > 0.0: continue
                        X += [1.0/float(x)]; X_arr[n] += [1.0/float(x)]
                        Y += [float(y)]; Y_arr[n] += [float(y)]
                    file.close()
                    if min(X) > x_min: x_min = min(X)
                
                X = X_arr[0]
                color_count = 0
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                sampArr = []
                for samp in range(1, int(max_n / 2) + 1):
                    print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                    Y, YErr = avgData(samp, Y_arr)
                    sampArr += [YErr]
                N_SAMP += [sampArr]
        for samp in range(len(N_SAMP[0])):
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,10))
            legend = []
            filename = "n_" + str(samp+1) + "_N"
            color_count = 0
            for N in range(len(N_Array)):
                print("plotting: %i sample: N: %i, J: %s, %i/%i          \r" % (samp+1, N_Array[N], str(J), J_cnt, len(valid_J)), end = "")
                YErr = N_SAMP[N][samp]
                # ploting
                if ABS:
                    subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                if REL: 
                    Xrel, Yrel, YErrrel = getRel(X, Y, YErr)
                    if min(Xrel) > x_min: x_min = min(Xrel)
                    subfig1.plot(Xrel, YErrrel, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
                legend += [Line2D([0], [0], label = "N = " + str(N_Array[N]), color = colors[color_count], ls = "solid", lw = line_width)]
                color_count += 1
                filename += "_" + str(N_Array[N])
                if color_count >= len(colors): break
            # saving
            handles, labels = plt.gca().get_legend_handles_labels()
            handles.extend(legend)
            subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_ylabel(r'$\sigma$ / $J_2$', fontsize = labelfontsize)
            vec = "einem Startvektor"
            if int(samp) > 1: vec = str(samp) + " Startvektoren"
            subfig1.set_title(r"$\sigma$" + " von $\\chi/N$ bei der QT mit Mittelung\nüber " + vec + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
            subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/sigma/X_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/sigma/X_sigma_" + abs_rel + "_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...                                        ")

def spin_gap_sigma():
    print("plotting sigma of spin gap ...")
    abs_rel = ""
    if ABS: abs_rel += "abs"
    if REL: abs_rel += "rel"
    for samp in range(1, int(max_n / 2) + 1):
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
        legend = []
        color_count = 0
        x_min = 0
        filename = "N"
        for N in N_Array:
            print("get data: samp: %i, N: %s                    \r" % (samp, N), end = "")
            N = str(N)
            X_arr = [[]] * max_n; Y_arr = [[]] * max_n; A_arr = [[]] * max_n
            for n in range(max_n):
                X_arr[n], Y_arr[n], A_arr[n] = getDataQT(N, n + 1)
            X, XErr = avgData(samp, X_arr); X = np.asarray(X)
            Y, YErr = avgData(samp, Y_arr); Y = np.asarray(Y); YErr = np.asarray(YErr)
            A, AErr = avgData(samp, A_arr); A = np.asarray(A); AErr = np.asarray(AErr)
            if J_SUM_SCALE: Y = Y / (1.0 + X); YErr = YErr / (1.0 + X)
            if ABS:
                subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig2.plot(X, AErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
            if REL: 
                Xrel, Yrel, YErrrel = getRel(X, Y, YErr)
                if min(Xrel) > x_min: x_min = min(Xrel)
                subfig1.plot(Xrel, YErrrel, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
                Xrel, Arel, AErrrel = getRel(X, A, AErr)
                if min(Xrel) > x_min: x_min = min(Xrel)
                subfig2.plot(Xrel, AErrrel, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
            legend += [Line2D([0], [0], label = "N = " + N, color = colors[color_count], ls = "solid", lw = line_width)]
            color_count += 1
            filename += "_" + str(N)
        # saving
        print("plotting data: samp: %i                    \r" % (samp), end = "")
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend(legend)
        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
        subfig2.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\sigma$ / $J_2$', fontsize = labelfontsize)
        subfig2.set_ylabel(r'$\sigma$ / $J_2$', fontsize = labelfontsize)
        vec = "einen Startvektor"
        if int(samp) > 1: vec = str(samp) + " Startvektoren"
        subfig1.set_title(r"$\sigma(\Delta_{SG})$" + " der Spinlückenenergie bei Mittelung über " + vec, fontsize = titlefontsize)
        subfig2.set_title(r"$\sigma(A)$" + " der Spinlückenenergie bei Mittelung über " + vec, fontsize = titlefontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        subfig2.tick_params(axis="both", which="major", labelsize=axisfontsize)
        subfig1.axhline(0, color = "grey")
        subfig2.axhline(0, color = "grey")
        subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig2.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig1.axvline(x=1.0, color='black', linestyle='--')
        subfig2.axvline(x=1.0, color='black', linestyle='--')
        if J_SUM_SCALE: filename = "J_SUM_SCALE_" + filename
        for end in [2.0, 2.5]:
            subfig1.set_xlim(0.0, end)
            subfig1.set_yscale('log')
            fig1.savefig("results/QT_stats/sigma/SG_k_sigma_" + abs_rel + "_samp_" + str(samp) + "_" + filename + "_" + str(start) + "_" + str(end) + ".png")
            fig1.savefig("results/QT_stats/sigma/SG_k_sigma_" + abs_rel + "_samp_" + str(samp) + "_" + filename + "_" + str(start) + "_" + str(end) + ".pdf")
            subfig2.set_xlim(0.0, end)
            subfig2.set_yscale('log')
            fig2.savefig("results/QT_stats/sigma/SG_A_sigma_" + abs_rel + "_samp_" + str(samp) + "_" + filename + "_" + str(start) + "_" + str(end) + ".png")
            fig2.savefig("results/QT_stats/sigma/SG_A_sigma_" + abs_rel + "_samp_" + str(samp) + "_" + filename + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.cla()
        plt.clf()
        plt.close(fig1)
        plt.close(fig2)
        plt.close('all')
        del fig1, subfig1
        del fig2, subfig2
        gc.collect()
    print("done...                              ")

if __name__ == "__main__":
    N_Array = [10, 12, 14, 16, 18, 20, 22]
    C_plot_n_for_each_N_sigma()
    print()
    # C_plot_N_for_each_n_sigma()
    # print()
    X_plot_n_for_each_N_sigma()
    print()
    # X_plot_N_for_each_n_sigma()
    # print()
    # spin_gap_sigma()
    # print()

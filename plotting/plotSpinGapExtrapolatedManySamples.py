import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
matplotlib.use('TkAgg')
import os
import sys
import numpy as np
import scipy.optimize
from typing import Tuple
import gc
import time
plt.rcParams['text.usetex'] = True

N_color = [("14", "red"), ("16", "blue"), ("18", "green"), ("20", "magenta"), ("22", "brown"), ("24", "purple"), ("26", "tomato")]
N_color_LOW = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color_HIGH = [("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta"), ("28", "brown"), ("30", "purple"), ("32", "tomato")]

peak_offset = 2000
fit_samples = 3 # 5
search_start_percent = 4/5
search_end_percent = 1/5
max_n = 30 # min = 1; max = 30

extrapolate_QT = True
extrapolate_ED = True

SHANK_ALG = False
EXP_FIT = False
ONE_OVER_N_FIT = True
EXPSILON_ALG = False

counter = 0

line_width = 1
marker_size = 4
alph = 0.1


def format_time(start_time: float, end_time: float) -> str:
    elapsed_time = end_time- start_time
    hours = int( elapsed_time / 3600 )
    minutes = int( (elapsed_time - hours * 3600) / 60 )
    seconds = elapsed_time - hours * 3600 - minutes * 60
    ret = ""
    if hours > 0: ret += str(hours) + " hours " + str(minutes) + " minutes "
    elif minutes > 0: ret += str(minutes) + " minutes "
    ret += str(seconds) + " seconds"
    return ret

def shank_Alg(A, n):
    Anp1 = A[n+1]; An = A[n]; Anm1 = A[n-1]
    return (Anp1 * Anm1 - An * An) / (Anp1 - 2 * An + Anm1) # Anp1 - ( (Anp1 - An) * (Anp1 - An) ) / ( (Anp1 - An) - (An - Anm1) )

def epsilon_Alg(A):
    length = len(A)
    eps = [[0] * (length + 1) for _ in range(length)]
    #for i in eps: print(i)
    for i in range(length):
        eps[i][0] = 0.0; eps[i][1] = A[i]
    for i in eps: print(i)
    for k in range(1, length + 1):
        for j in range(length - k):
            #print("j: %i, k: %i" % (j, k))
            eps[j][k+1] = eps[j+1][k-1] + 1.0 / (eps[j+1][k] - eps[j][k])
            #for i in eps: print(i)
    for i in eps: print(i)
    epslim = -1.0
    if (length+1) % 2 == 0: epslim = eps[0][length]
    else: epslim = eps[1][length-1]
    print(epslim)
    return epslim

def expFunc_offset(x: float, A: float, k: float, x_0: float) -> float:
    return A * np.exp(k * x) + x_0

def one_over_N(x: float, A: float, b: float) -> float:
    return A / x + b

def expFunc(x: float, A: float, k: float) -> float:
    return x * A * np.exp(k * x)

def extrapolation_alg_raw(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]

        # fitting exp to data
        if EXP_FIT:
            X = [x for x in range(1, len(A)+1)]
            fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
            subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
            counter += 1
            try:
                X = np.array(X); A = np.array(A)
                params, cv = scipy.optimize.curve_fit(expFunc_offset, X, A, (1.0, 0.1, 0.1))
                A_param, k_param, x_0_param = params
                X = np.linspace(0, X[-1], 100)
                Y = expFunc_offset(X, A_param, k_param, x_0_param)
                subfig3.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = "black")
                outdata += [x_0_param]
            except RuntimeError:
                print("fit failed")
                outdata += [-1]
            fig3.savefig("results/" + "extrapolation_temp/data_" + str(counter) + "_.png")
            plt.close(fig3)

        # fitting 1/N to data
        elif ONE_OVER_N_FIT:
            X = [x for x in range(1, len(A)+1)]
            fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
            counter += 1
            try:
                X = np.array(X); A = np.array(A)
                X = 1 / X
                subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
                # params, cv = scipy.optimize.curve_fit(one_over_N, X, A, (0.1, 0.1))
                # A_param, b_param = params
                A_param, b_param = np.polyfit(X, A, 1)
                X = np.linspace(0.0, X[0], 100)
                # Y = one_over_N(X, A_param, b_param)
                Y = A_param * X + b_param
                subfig3.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = "black")
                outdata += [abs(b_param)]
            except RuntimeError:
                print("fit failed")
                outdata += [-1]
            fig3.savefig("results/" + "extrapolation_temp/data_" + str(counter) + "_.png")
            plt.close(fig3)

        # shank algorithm
        elif SHANK_ALG:
            S = []
            for j in range(1, len(A)-1): S += [shank_Alg(A, j)]
            outdata += [S[-1]]

        # epsilon algorithm
        elif EXPSILON_ALG:
            outdata += [epsilon_Alg(A)]

    return outdata

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
    file = open("results/" + N + "/data/spin_gap_data/" + str(n+1) + "/" + filename, 'r')
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

def getDataQT(N: int, n: int):
    X = []; Y = [];  A = []
    for filename in os.listdir("results/" + N + "/data/spin_gap_data/" + str(n) + "/"):
        if filename[len(filename)-6:] != "QT.txt": continue
        J = filename[len("X_J"):-len("QT.txt")]
        A, k = get_spin_gap(n, N, J, filename)
        X += [float(J)]
        Y += [float(k)]
        A += [float(A)]
    X, Y, A = sort_data(X, Y, A)
    return X, Y, A

def getDataED(N: int, n: int):
    X = []; Y = [];  A = []
    for filename in os.listdir("results/" + N + "/data/spin_gap_data/" + str(n) + "/"):
        if filename[len(filename)-6:] != "ED.txt": continue
        J = filename[len("X_J"):-len("ED.txt")]
        A, k = get_spin_gap(n, N, J, filename)
        X += [float(J)]
        Y += [float(k)]
        A += [0]
    X, Y, A = sort_data(X, Y, A)
    return X, Y, A

def getExtrapolatedDataQT(N_color):
    print("getting extrapolating data (QT)")
    extrapolated_data = [] * max_n
    for n in range(1, max_n + 1):
        print("n = %i / %i\r" %(n, max_n), end="")
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        extrapolationArray = []
        used_N = "N"
        for N, C in N_color:
            try:
                X, Y, A = getDataQT(int(N), n)
                extrapolationArray += [Y]
                used_N += "_" + N
                subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 4, marker = "o", color = C, label = "N = " + N)
            except:
                print("QT data for N = %s no available" % N)

        extrapolated_data[n-1] = extrapolation_alg_raw(extrapolationArray)

        subfig1.plot(X, extrapolated_data[n-1], lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = "black", label = "Extrapolation")

        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 25)
        subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$ mit einem Startvektor', fontsize = 25)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

        fig1.savefig("results/" + "SE_raw_extrap_QT_" + used_N + ".png")
        plt.close(fig1)
    return extrapolated_data

def getExtrapolatedDataED(N_color):
    print("getting extrapolating data (ED)")
    extrapolated_data = []
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    extrapolationArray = []
    used_N = "N"
    for N, C in N_color:
        try:
            X, Y, A = getDataED(int(N), 1)
            extrapolationArray += [Y]
            used_N += "_" + N
            subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 4, marker = "o", color = C, label = "N = " + N)
        except:
            print("ED data for N = %s no available" % N)

    extrapolated_data = extrapolation_alg_raw(extrapolationArray)

    subfig1.plot(X, extrapolated_data, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = "black", label = "Extrapolation")

    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
    subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 25)
    subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$ mit einem Startvektor', fontsize = 25)

    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

    fig1.savefig("results/" + "SE_raw_extrap_QT_" + used_N + ".png")
    plt.close(fig1)
    return extrapolated_data

def avgData(samp, inputdata):
    data_to_avg = []
    for i in range(0, max_n - samp, samp):
        data_raw = []
        for j in range(i, i + samp):
            data_raw += [inputdata[j]]
        Y_data = []
        for k in range(0, len(data_raw[0])):
            Y_data_raw = []
            for Yd in data_raw:
                Y_data_raw += [Yd[k]]
            Y_data_raw = np.array(Y_data_raw)
            Y_data += [Y_data_raw.mean()]
        data_to_avg += [Y_data]
    Y = []; YErr = []
    for k in range(0, len(data_to_avg[0])):
        Y_data_raw = []
        for Yd in data_to_avg:
            Y_data_raw += [Yd[k]]
        Y_data_raw = np.array(Y_data_raw)
        Y += [Y_data_raw.mean()]; YErr +=[Y_data_raw.std()]
    return Y, YErr

def plotExtrapolatedData(N_color):
    print("plotting extrapolating data")
    extrapolated_data_QT = getExtrapolatedDataQT(N_color)
    extrapolated_data_ED = getExtrapolatedDataED(N_color)
    for samp in range(1, max_n / 2 + 1):
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        used_N_QT = "N"; used_N_ED = "N"; N_arr = [0] * 35

        # ED #
        for N, C in N_color:        
            try:
                X_arr = [] * max_n; Y_arr = [] * max_n; A_arr = [] * max_n
                for n in range(max_n):
                    X_arr[n], Y_arr[n], A_arr[n] = getDataED(int(N), n + 1)
                X, XErr = avgData(samp, X_arr); X = np.asarray(X)
                Y, YErr = avgData(samp, Y_arr); Y = np.asarray(Y); YErr = np.asarray(YErr)
                subfig1.plot(X, Y, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = C)#, label = "N = " + N)
                subfig1.fill_between(X, Y - YErr, Y + YErr, color = C, alpha = alph)
                used_N_ED += "_" + N; N_arr[int(N)] = 1
            except:
                print("cannot plot ED data for N = %s" % N)
        Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_ED)
        Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
        subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = "black")#, label = "Extrapolation")
        subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = "black", alpha = alph)

        # QT #
        for N, C in N_color:        
            try:
                X_arr = [] * max_n; Y_arr = [] * max_n; A_arr = [] * max_n
                for n in range(max_n):
                    X_arr[n], Y_arr[n], A_arr[n] = getDataQT(int(N), n + 1)
                X, XErr = avgData(samp, X_arr); X = np.asarray(X)
                Y, YErr = avgData(samp, Y_arr); Y = np.asarray(Y); YErr = np.asarray(YErr)
                subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = C)#, label = "N = " + N)
                subfig1.fill_between(X, Y - YErr, Y + YErr, color = C, alpha = alph)
                used_N_QT += "_" + N; N_arr[int(N)] = 1
            except:
                print("cannot plot QT data for N = %s" % N)

        Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT)
        Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
        subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = "black")#, label = "Extrapolation")
        subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = "black", alpha = alph)

        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 25)
        if samp == 1: vec_string = "einem Startvektor"
        else: vec_string = str(samp) + " Startvektoren"
        subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$' + "\nmit " + vec_string + " und je " + str(int(max_n/samp)) + " Mittlungen", fontsize = 25)

        legend = []
        for i in range(len(N_arr)):
            N_index = N_arr[i]
            if N_index == 1:
                N = str(i)
                C = ""
                for n, c in N_color:
                    if n == N: C = c
                legend += [Line2D([0], [0], label = "N = " + N, color = C)]
        legend += [Line2D([0], [0], label = "Extrapolation", color = "black")]
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend(legend)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = 20)
        subfig1.tick_params(axis="both", which="major", labelsize=25)

        fig1.savefig("results/" + "SE_extrap_ED_" + used_N_ED + "_QT_" + used_N_QT  + "_max_" + str(max_n) + "_samp_" + str(samp) + ".png")
        fig1.savefig("results/" + "SE_extrap_ED_" + used_N_ED + "_QT_" + used_N_QT  + "_max_" + str(max_n) + "_samp_" + str(samp) + ".pdf")
        plt.close(fig1)

if __name__ == "__main__":
    start_time = time.time()

    # if len(sys.argv) > 1:
    #     regime = sys.argv[1]
    #     if regime == "low": N_color = N_color_LOW#; print("low regime")
    #     elif regime == "high": N_color = N_color_HIGH; no_ED = True#; print("high regime")
    #     else: N_color = N_color_LOW; print("default low (wrong args)")
    # else: print("default (no args)")

    N_color = [("12", "red"), ("14", "blue"), ("16", "green"), ("18", "magenta"), ("20", "brown"), ("22", "purple"), ("24", "tomato")]

    plotExtrapolatedData(N_color)
    print()

    end_time = time.time()

    print("done; this took %s" % format_time(start_time, end_time))

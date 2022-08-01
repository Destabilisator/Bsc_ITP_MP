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

N_color = [("12", "red"), ("14", "blue"), ("16", "green"), ("18", "magenta"), ("20", "brown"), ("22", "purple"), ("24", "tomato")]
N_color = [("10", "red"), ("12", "blue"), ("14", "green"), ("16", "magenta"), ("18", "brown"), ("20", "purple"), ("22", "tomato")]
N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple")]#, ("18", "tomato")]

colors = ["red", "blue", "green", "magenta", "brown", "purple", "tomato"]

peak_offset = 2000
fit_samples = 3 # 5
search_start_percent = 4/5
search_end_percent = 1/5
max_n = 5 # min = 1; max = 30

extrapolate_QT = True
extrapolate_ED = True

titlefontsize = 35
labelfontsize = 30
legendfontsize = 30
axisfontsize = 25

# in final output
EXP_FIT = False
ONE_OVER_N_FIT = True
ONE_OVER_N2_FIT = False
ONE_OVER_SQRTN_FIT = False
SHANK_ALG = False
EXPSILON_ALG = False

# in multiple extrapolations output
SHOW_EXP_FIT = False
SHOW_ONE_OVER_N_FIT = True
SHOW_ONE_OVER_N2_FIT = True
SHOW_ONE_OVER_SQRTN_FIT = True
SHOW_SHANK_ALG = True
SHOW_EXPSILON_ALG = True

MAKE_NEW = False

counter = 0

line_width = 4
marker_size = 5
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
    #for i in eps: print(i)
    for k in range(1, length + 1):
        for j in range(length - k):
            #print("j: %i, k: %i" % (j, k))
            eps[j][k+1] = eps[j+1][k-1] + 1.0 / (eps[j+1][k] - eps[j][k])
            #for i in eps: print(i)
    #for i in eps: print(i)
    epslim = -1.0
    if (length+1) % 2 == 0: epslim = eps[0][length]
    else: epslim = eps[1][length-1]
    #print(epslim)
    return epslim

def expFunc_offset(x: float, A: float, k: float, x_0: float) -> float:
    return A * np.exp(k * x) + x_0

def one_over_N(x: float, A: float, b: float) -> float:
    return A / x + b

def expFunc(x: float, A: float, k: float) -> float:
    return x * A * np.exp(k * x)

def exp_fit_extrap(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]
        # fitting exp to data
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
    return outdata

def sqrtN_fit_extrap(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]
        # fitting 1/sqrt(N) to data
        X = [x for x in range(1, len(A)+1)]
        fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
        counter += 1
        try:
            X = np.array(X); A = np.array(A)
            X = 1 / np.sqrt(X)
            subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
            # params, cv = scipy.optimize.curve_fit(one_over_N, X, A, (0.1, 0.1))
            # A_param, b_param = params
            A_param, b_param = np.polyfit(X, A, 1)# , w = 1/X + 1/2 * 1/X**2  + 1/6 * 1/X**3
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
    return outdata

def N_fit_extrap(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]
        # fitting 1/N to data
        X = [x for x in range(1, len(A)+1)]
        fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
        counter += 1
        try:
            X = np.array(X); A = np.array(A)
            X = 1 / X
            subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
            # params, cv = scipy.optimize.curve_fit(one_over_N, X, A, (0.1, 0.1))
            # A_param, b_param = params
            A_param, b_param = np.polyfit(X, A, 1)# , w = 1/X + 1/2 * 1/X**2  + 1/6 * 1/X**3
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
    return outdata

def N2_fit_extrap(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]
        # fitting 1/N**2 to data
        X = [x for x in range(1, len(A)+1)]
        fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
        counter += 1
        try:
            X = np.array(X); A = np.array(A)
            X = 1 / X**2
            subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
            # params, cv = scipy.optimize.curve_fit(one_over_N, X, A, (0.1, 0.1))
            # A_param, b_param = params
            A_param, b_param = np.polyfit(X, A, 1)# , w = 1/X + 1/2 * 1/X**2  + 1/6 * 1/X**3
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
    return outdata

def expsilon_extrap(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]
        # epsilon algorithm
        outdata += [epsilon_Alg(A)]
    return outdata

def shank_extrap(A_raw):
    global counter
    len_N = len(A_raw); len_Y = len(A_raw[0])
    outdata = []
    for i in range(len_Y):
        A = []
        for j in range(len_N): A += [A_raw[j][i]]
        # shank algorithm
        S = []
        for j in range(1, len(A)-1): S += [shank_Alg(A, j)]
        outdata += [S[-1]]
    return outdata

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
                A_param, b_param = np.polyfit(X, A, 1)# , w = 1/X + 1/2 * 1/X**2  + 1/6 * 1/X**3
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
        
        # fitting 1/N**2 to data
        elif ONE_OVER_N2_FIT:
            X = [x for x in range(1, len(A)+1)]
            fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
            counter += 1
            try:
                X = np.array(X); A = np.array(A)
                X = 1 / X**2
                subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
                # params, cv = scipy.optimize.curve_fit(one_over_N, X, A, (0.1, 0.1))
                # A_param, b_param = params
                A_param, b_param = np.polyfit(X, A, 1)# , w = 1/X + 1/2 * 1/X**2  + 1/6 * 1/X**3
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

        # fitting 1/sqrt(n) to data
        elif ONE_OVER_SQRTN_FIT:
            X = [x for x in range(1, len(A)+1)]
            fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
            counter += 1
            try:
                X = np.array(X); A = np.array(A)
                X = 1 / np.sqrt(X)
                subfig3.plot(X, A, lw = 0, ls = "dashed", markersize = 5, marker = "o", color = "black")
                # params, cv = scipy.optimize.curve_fit(one_over_N, X, A, (0.1, 0.1))
                # A_param, b_param = params
                A_param, b_param = np.polyfit(X, A, 1)# , w = 1/X + 1/2 * 1/X**2  + 1/6 * 1/X**3
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
    return X, Y, A

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

def getExtrapolatedDataQT(N_color):
    print("getting extrapolated data (QT)")
    extrapolated_data = [[]] * max_n
    for n in range(1, max_n + 1):
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        extrapolationArray = []
        used_N = "N"
        for N, C in N_color:
            print("n = %i / %i & N = %s (get data)          \r" %(n, max_n, N), end="")
            try:
                X, Y, A = getDataQT(N, n)
                extrapolationArray += [Y]
                used_N += "_" + N
                subfig1.plot(X, Y, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = C, label = "N = " + N)
            except:
                print("QT data for N = %s no data available" % N)

        print("n = %i / %i (extrapolating data)          \r" %(n, max_n), end="")
        extrapolated_data[n-1] = extrapolation_alg_raw(extrapolationArray)

        subfig1.plot(X, extrapolated_data[n-1], lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = "black", label = "Extrapolation")

        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = labelfontsize)
        subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$ mit einem Startvektor', fontsize = titlefontsize)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

        fig1.savefig("results/" + "SG/SG_raw_extrap_QT_" + used_N + "_" + str(n) + ".png")
        plt.close(fig1)
    return extrapolated_data

def newEDdata(N: str):
    X = []; Y = [];  A = []
    for filename in os.listdir("results/" + N + "/data/spin_gap_data/1/"):
        if filename[len(filename)-6:] != "ED.txt": continue
        J = filename[len("X_J"):-len("ED.txt")]
        amp, k = get_spin_gap(1, N, J, filename)
        X += [float(J)]
        Y += [float(k)]
        A += [float(amp)]
    X, Y, A = sort_data(X, Y, A)
    file = open("results/" + N + "/data/spin_gap_data/ED.txt", "w")
    for i in range(len(X)):
        file.write(f"{X[i]}\t{Y[i]}\t{A[i]}\n")
    file.close()
    return X, Y, A

def getDataED(N: str):
    if MAKE_NEW:
        X = []; Y = [];  A = []
        X, Y, A = newEDdata(N)
    else:
        try:
            file = open("results/" + N + "/data/spin_gap_data/ED.txt", 'r')
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
            X, Y, A = newEDdata(N)
    return X, Y, A

def getExtrapolatedDataED(N_color):
    print("getting extrapolated data (ED)")
    extrapolated_data = []
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    extrapolationArray = []
    used_N = "N"
    for N, C in N_color:
        print("N = %s (get data)          \r" %(N), end="")
        try:
            X, Y, A = getDataED(N)
            extrapolationArray += [Y]
            used_N += "_" + N
            subfig1.plot(X, Y, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = C, label = "N = " + N)
        except:
            print("ED data for N = %s no data available" % N)

    print("ED (extrapolating data)          \r", end="")
    extrapolated_data = extrapolation_alg_raw(extrapolationArray)

    subfig1.plot(X, extrapolated_data, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = "black", label = "Extrapolation")

    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = labelfontsize)
    subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$ mit einem Startvektor', fontsize = titlefontsize)

    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

    fig1.savefig("results/" + "SG/SG_raw_extrap_ED_" + used_N + ".png")
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
    print("getting extrapolated data")
    extrapolated_data_QT = getExtrapolatedDataQT(N_color)
    extrapolated_data_ED = getExtrapolatedDataED(N_color)
    print("plotting extrapolated data")
    for samp in range(1, int(max_n / 2) + 1):
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        used_N_QT = "N"; used_N_ED = "N"; N_arr = [0] * 35

        # ED #
        for N, C in N_color:        
            try:
                X, Y, A = getDataED(N)
                X = np.asarray(X); Y = np.asarray(Y)
                subfig1.plot(X, Y, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = C, alpha = 0.5)#, label = "N = " + N)
                used_N_ED += "_" + N; N_arr[int(N)] = 1
            except:
                print("cannot plot ED data for N = %s" % N)
        subfig1.plot(X, extrapolated_data_ED, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = "black", alpha = 0.75)#, label = "Extrapolation")
        
        # QT #
        for N, C in N_color:        
            try:
                X_arr = [[]] * max_n; Y_arr = [[]] * max_n; A_arr = [[]] * max_n
                for n in range(max_n):
                    X_arr[n], Y_arr[n], A_arr[n] = getDataQT(N, n + 1)
                X, XErr = avgData(samp, X_arr); X = np.asarray(X)
                Y, YErr = avgData(samp, Y_arr); Y = np.asarray(Y); YErr = np.asarray(YErr)
                subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = C, alpha = 0.5)#, label = "N = " + N)
                subfig1.fill_between(X, Y - YErr, Y + YErr, color = C, alpha = alph)
                used_N_QT += "_" + N; N_arr[int(N)] = 1
            except:
                print("cannot plot QT data for N = %s" % N)
        Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT)
        Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
        subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = "black", alpha = 0.75)#, label = "Extrapolation")
        subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = "black", alpha = alph)

        # saving plots #
        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\Delta E_{gap}$ / $J_2$', fontsize = labelfontsize)
        if samp == 1: vec_string = "einen Startvektor"
        else: vec_string = str(samp) + " Startvektoren"
        subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$' + "\nmit " + str(int(max_n/samp)) + " Mittlungen über " + vec_string, fontsize = titlefontsize)

        legend = []
        for i in range(len(N_arr)):
            N_index = N_arr[i]
            if N_index == 1:
                N = str(i)
                C = ""
                for n, c in N_color:
                    if n == N: C = c
                legend += [Line2D([0], [0], label = "N = " + N, color = C, ls = "solid", lw = line_width)]
        legend += [Line2D([0], [0], label = "Extrapolation", color = "black", ls = "solid", lw = line_width)]
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend(legend)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

        fig1.savefig("results/" + "SG/SG_extrap_ED_" + used_N_ED + "_QT_" + used_N_QT  + "_max_" + str(max_n) + "_samp_" + str(samp) + ".png")
        fig1.savefig("results/" + "SG/SG_extrap_ED_" + used_N_ED + "_QT_" + used_N_QT  + "_max_" + str(max_n) + "_samp_" + str(samp) + ".pdf")
        plt.close(fig1)

def different_extrapolations(N_color):
    print("extrapolating with different algorithms")
    for samp in range(1, int(max_n / 2) + 1):
        if samp != 1: continue
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        used_N_QT = "N"; used_N_ED = "N"; N_arr = [0] * 35

        # ED #
        print("ED getting data          ")
        extrapolationArray = []
        for N, C in N_color:        
            try:
                X, Y, A = getDataED(N)
                extrapolationArray += [Y]
                X = np.asarray(X); Y = np.asarray(Y)
                subfig1.plot(X, Y, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = "black", alpha = 0.5)#, label = "N = " + N)
                used_N_ED += "_" + N; N_arr[int(N)] = 1
            except:
                print("cannot plot ED data for N = %s" % N)

        color_count = 0
        if SHOW_EXP_FIT:
            print("ED extrapolating with exp          \r", end = "")
            extrapolated_data_ED_EXP_FIT = exp_fit_extrap(extrapolationArray)
            subfig1.plot(X, extrapolated_data_ED_EXP_FIT, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75); color_count += 1
        if SHOW_ONE_OVER_SQRTN_FIT:
            print("ED extrapolating with sqrt          \r", end = "")
            extrapolated_data_ED_ONE_OVER_SQRTN_FIT = sqrtN_fit_extrap(extrapolationArray)
            subfig1.plot(X, extrapolated_data_ED_ONE_OVER_SQRTN_FIT, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75); color_count += 1
        if SHOW_ONE_OVER_N_FIT:
            print("ED extrapolating with N          \r", end = "")
            extrapolated_data_ED_ONE_OVER_N_FIT = N_fit_extrap(extrapolationArray)
            subfig1.plot(X, extrapolated_data_ED_ONE_OVER_N_FIT, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75); color_count += 1
        if SHOW_ONE_OVER_N2_FIT:
            print("ED extrapolating with N**2          \r", end = "")
            extrapolated_data_ED_ONE_OVER_N2_FIT = N2_fit_extrap(extrapolationArray)
            subfig1.plot(X, extrapolated_data_ED_ONE_OVER_N2_FIT, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75); color_count += 1
        if SHOW_SHANK_ALG:
            print("ED extrapolating with shank          \r", end = "")
            extrapolated_data_ED_SHANK_ALG = shank_extrap(extrapolationArray)
            subfig1.plot(X, extrapolated_data_ED_SHANK_ALG, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75); color_count += 1
        if SHOW_EXPSILON_ALG:
            print("ED extrapolating with epsilon          \r", end = "")
            extrapolated_data_ED_EPSILON_ALG = expsilon_extrap(extrapolationArray)
            subfig1.plot(X, extrapolated_data_ED_EPSILON_ALG, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75); color_count += 1
        print("ED extrapolation done          \r", end = "")

        # QT #
        print("QT getting data          ")
        for N, C in N_color:        
            try:
                X_arr = [[]] * max_n; Y_arr = [[]] * max_n; A_arr = [[]] * max_n
                for n in range(max_n):
                    X_arr[n], Y_arr[n], A_arr[n] = getDataQT(N, n + 1)
                X, XErr = avgData(samp, X_arr); X = np.asarray(X)
                Y, YErr = avgData(samp, Y_arr); Y = np.asarray(Y); YErr = np.asarray(YErr)
                subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = "black", alpha = 0.5)#, label = "N = " + N)
                subfig1.fill_between(X, Y - YErr, Y + YErr, color = "black", alpha = alph)
                used_N_QT += "_" + N; N_arr[int(N)] = 1
            except:
                print("cannot plot QT data for N = %s          " % N)

        extrapolated_data_QT_EXP = []
        extrapolated_data_QT_ONE_OVER_SQRTN_FIT = []
        extrapolated_data_QT_OVER_N_FIT = []
        extrapolated_data_QT_OVER_N2_FIT = []
        extrapolated_data_QT_SHANK_ALG = []
        extrapolated_data_QT_EPSILON_ALG = []

        for n in range(1, max_n + 1):
            extrapolationArray = []
            used_N = "N"
            for N, C in N_color:
                print("n = %i / %i & N = %s (get data)          \r" %(n, max_n, N), end="")
                try:
                    X, Y, A = getDataQT(N, n)
                    extrapolationArray += [Y]
                    used_N += "_" + N
                except:
                    print("QT data for N = %s no data available          " % N)

            print("n = %i / %i (extrapolating data)          \r" %(n, max_n), end="")
            
            if SHOW_EXP_FIT:
                print("n = %i / %i extrapolating with exp          \r", end = "")
                extrapolated_data_QT_EXP += [exp_fit_extrap(extrapolationArray)]
            if SHOW_ONE_OVER_SQRTN_FIT:
                print("n = %i / %i extrapolating with sqrt          \r" %(n, max_n), end = "")
            if SHOW_ONE_OVER_N_FIT:
                extrapolated_data_QT_ONE_OVER_SQRTN_FIT += [sqrtN_fit_extrap(extrapolationArray)]
                print("n = %i / %i extrapolating with N          \r" %(n, max_n), end = "")
                extrapolated_data_QT_OVER_N_FIT += [N_fit_extrap(extrapolationArray)]
            if SHOW_ONE_OVER_N2_FIT:
                print("n = %i / %i extrapolating with N**2          \r" %(n, max_n), end = "")
                extrapolated_data_QT_OVER_N2_FIT += [N2_fit_extrap(extrapolationArray)]
            if SHOW_SHANK_ALG:
                print("n = %i / %i extrapolating with shank          \r" %(n, max_n), end = "")
                extrapolated_data_QT_SHANK_ALG += [shank_extrap(extrapolationArray)]
            if SHOW_EXPSILON_ALG:
                print("n = %i / %i extrapolating with epsilon          \r" %(n, max_n), end = "")
                extrapolated_data_QT_EPSILON_ALG += [expsilon_extrap(extrapolationArray)]
        print("QT extrapolation done          \r", end = "")

        color_count = 0
        if SHOW_EXP_FIT:
            print("plotting QT extrapolating with exp          \r", end = "")
            Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT_EXP)
            Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
            subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75)
            subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = colors[color_count], alpha = alph); color_count += 1
        if SHOW_ONE_OVER_SQRTN_FIT:
            print("plotting QT QT extrapolating with sqrt          \r", end = "")
            Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT_ONE_OVER_SQRTN_FIT)
            Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
            subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75)
            subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = colors[color_count], alpha = alph); color_count += 1
        if SHOW_ONE_OVER_N_FIT:
            print("plotting QT QT extrapolating with N          \r", end = "")
            Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT_OVER_N_FIT)
            Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
            subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75)
            subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = colors[color_count], alpha = alph); color_count += 1
        if SHOW_ONE_OVER_N2_FIT:
            print("plotting QT QT extrapolating with N**2          \r", end = "")
            Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT_OVER_N2_FIT)
            Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
            subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75)
            subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = colors[color_count], alpha = alph); color_count += 1
        if SHOW_SHANK_ALG:
            print("plotting QT QT extrapolating with shank          \r", end = "")
            Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT_SHANK_ALG)
            Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
            subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75)
            subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = colors[color_count], alpha = alph); color_count += 1
        if SHOW_EXPSILON_ALG:
            print("plotting QT QT extrapolating with epsilon          \r", end = "")
            Y_etrap, YErr_etrap = avgData(samp, extrapolated_data_QT_EPSILON_ALG)
            Y_etrap = np.asarray(Y_etrap); YErr_etrap = np.asarray(YErr_etrap)
            subfig1.plot(X, Y_etrap, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = colors[color_count], alpha = 0.75)
            subfig1.fill_between(X, Y_etrap - YErr_etrap, Y_etrap + YErr_etrap, color = colors[color_count], alpha = alph); color_count += 1

        # saving plots #
        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\Delta E_{gap}$ / $J_2$', fontsize = labelfontsize)
        if samp == 1: vec_string = "einen Startvektor"
        else: vec_string = str(samp) + " Startvektoren"
        subfig1.set_title(r'Spinlücken-Energien $\Delta E_{gap}$' + "\nmit " + str(int(max_n/samp)) + " Mittlungen über " + vec_string, fontsize = titlefontsize)

        legend = []
        color_count = 0
        legend += [Line2D([0], [0], label = "Daten", color = "black", ls = "solid", lw = line_width, alpha = 0.5)]
        if SHOW_EXP_FIT: legend += [Line2D([0], [0], label = "exp-Fit", color = colors[color_count], ls = "solid", lw = line_width)]; color_count += 1
        if SHOW_ONE_OVER_SQRTN_FIT: legend += [Line2D([0], [0], label = r"1/$\sqrt{N}$-Fit", color = colors[color_count], ls = "solid", lw = line_width)]; color_count += 1
        if SHOW_ONE_OVER_N_FIT: legend += [Line2D([0], [0], label = r"1/$N$-Fit", color = colors[color_count], ls = "solid", lw = line_width)]; color_count += 1
        if SHOW_ONE_OVER_N2_FIT: legend += [Line2D([0], [0], label = r"1/$N^2$-Fit", color = colors[color_count], ls = "solid", lw = line_width)]; color_count += 1
        if SHOW_SHANK_ALG: legend += [Line2D([0], [0], label = "Shanks-Tranformation", color = colors[color_count], ls = "solid", lw = line_width)]; color_count += 1
        if SHOW_EXPSILON_ALG: legend += [Line2D([0], [0], label = r"$\varepsilon$-Algorithmus", color = colors[color_count], ls = "solid", lw = line_width)]; color_count += 1
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend(legend)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

        fig1.savefig("results/" + "SG/SG_extrap_methods_ED_" + used_N_ED + "_QT_" + used_N_QT  + "_max_" + str(max_n) + "_samp_" + str(samp) + ".png")
        fig1.savefig("results/" + "SG/SG_extrap_methods_ED_" + used_N_ED + "_QT_" + used_N_QT  + "_max_" + str(max_n) + "_samp_" + str(samp) + ".pdf")
        plt.close(fig1)

if __name__ == "__main__":
    start_time = time.time()

    # if len(sys.argv) > 1:
    #     regime = sys.argv[1]
    #     if regime == "low": N_color = N_color_LOW#; print("low regime")
    #     elif regime == "high": N_color = N_color_HIGH; no_ED = True#; print("high regime")
    #     else: N_color = N_color_LOW; print("default low (wrong args)")
    # else: print("default (no args)")

    # plotExtrapolatedData(N_color)
    # print()

    N_color = [("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple")]#, ("18", "tomato")]

    # plotExtrapolatedData(N_color)
    # print()
    
    different_extrapolations(N_color)
    print()

    end_time = time.time()

    print("done; this took %s" % format_time(start_time, end_time))

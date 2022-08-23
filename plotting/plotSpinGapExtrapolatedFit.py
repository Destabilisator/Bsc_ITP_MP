from tkinter import PIESLICE
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
N_color_FULL = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato"),
                ("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta"), ("28", "brown"), ("30", "purple"), ("32", "tomato")]

N_color = [("12", "red"), ("14", "blue"), ("16", "green"), ("18", "magenta"), ("20", "brown"), ("22", "purple"), ("24", "tomato")]
N_color = [("10", "red"), ("12", "blue"), ("14", "green"), ("16", "magenta"), ("18", "brown"), ("20", "purple"), ("22", "tomato")]
# N_color = [("12", "blue"), ("14", "green"), ("16", "magenta"), ("18", "brown"), ("20", "purple"), ("22", "tomato")]

N_vals_QT = ["6", "8", "10", "12", "14", "16", "18", "20", "22"]
N_vals_ED = ["6", "8", "10", "12", "14", "16"]

colors = ["red", "blue", "green", "magenta", "brown", "purple", "tomato"]

peak_offset = 2000
fit_samples = 3 # 5
search_start_percent = 4/5
search_end_percent = 1/5
max_n = 30 # min = 1; max = 30
epsilon = 10**(-7)

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

# force generate new data
MAKE_NEW = False

# scale data with J1 + J2
J_SUM_SCALE = True

# change output
LINEAR_FIT = True
ONE_OVER_N_FIT = True
ONE_OVER_N2_FIT = True
ONE_OVER_SQRTN_FIT = True

# fehlerbalken der QT ausblenden
QT_WITH_ERROR = False


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
    return x * A * np.exp(k * x)

def linFunc(x: float, m: float, b: float) -> float:
    return m * x + b

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

def get_spin_gap(n: int, N: str, J: str, filename: str) -> Tuple[float, float]:
    file = open("results/" + N + "/data/spin_gap_data/" + str(n) + "/" + filename, 'r')
    ED_QT = filename[len(filename)-6:-4]
    lines = file.readlines()
    X = []; Y = []
    for i in range(5,len(lines)):
        arr = lines[i].split("\t")
        x = arr[0]
        y = arr[1]
        if float(y) < epsilon: continue
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
    # X = np.asarray(X)
    if len(X) == 0: raise ValueError('n = %i: QT data does not exist' % n)
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

def generateSGDataQT(N_color):
    print("generating SG data (QT)          ")
    for n in range(1, max_n + 1):
        for N, C in N_color:
            print("n = %i / %i & N = %s (generate data)          \r" %(n, max_n, N), end="")
            try:
                X, Y, A = getDataQT(N, n)
            except:
                print("QT data for N = %s and n = %i not available" %(N, n))
    print("done...                              ")

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
    if len(X) == 0: raise ValueError('ED data does not exist')
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

def generateSGDataED(N_color):
    print("generating SG data (ED)          ")
    for N, C in N_color:
        print("N = %s (generate data)          \r" %(N), end="")
        try:
            X, Y, A = getDataED(N)
        except:
            print("ED data for N = %s not available" % N)
    print("done...                              ")

def getSGDataQT():
    out = []
    outErr = []
    for N in N_vals_QT:
        X_out = []; Y_out = []; A_out = []
        # X_out = np.asarray(X_out); Y_out = np.asarray(Y_out); A_out = np.asarray(A_out)
        for n in range(1, max_n + 1):
            X, Y, A = getDataQT(N, n)
            X = np.asarray(X); Y = np.asarray(Y); A = np.asarray(A)
            if J_SUM_SCALE: Y = Y / (1.0 + X)
            X_out += [X]; Y_out += [Y]; A_out += [A]
        X = []; Y = []; YErr = []; A = []; AErr = []
        for pos in range(len(X_out[0])):
                y = []; a = []
                for n in range(max_n):
                    if len(X_out[n]) == 0: continue
                    y += [Y_out[n][pos]]
                    a += [A_out[n][pos]]
                y = np.array(y); a = np.array(a)
                X += [X_out[0][pos]]; Y += [y.mean()]; YErr += [y.std()]; A += [a.mean()]; AErr += [a.std()/a.mean()]
        out += [Y]
        outErr += [YErr]
    return out, outErr

def getSGDataED():
    out = []
    for N in N_vals_ED:
        X, Y, A = getDataED(N)
        X = np.asarray(X); Y = np.asarray(Y); A = np.asarray(A)
        if J_SUM_SCALE: Y = Y / (1.0 + X)
        out += [Y]
    return out

if __name__ == "__main__":
    start_time = time.time()

    if MAKE_NEW: generateSGDataQT(N_color_FULL); generateSGDataED(N_color_FULL)

    print("getting data...")
    ED_data = getSGDataED(); QT_data, QT_data_Err = getSGDataQT()
    X, Y, A = getDataQT("6", 1)

    for J in range(len(X)):
        J_val = str(X[J])
        print("J = %s          \r" %J_val, end="")

        figData_ED_X = []; figData_QT_X = []
        figData_ED_Y = []; figData_QT_Y = []; figData_QT_YErr = []
        for n in range(len(N_vals_QT)):
            figData_QT_X += [int(N_vals_QT[n])]
            figData_QT_Y += [QT_data[n][J]]
            figData_QT_YErr += [QT_data_Err[n][J]*100]
        for n in range(len(N_vals_ED)):
            figData_ED_X += [int(N_vals_ED[n])]
            figData_ED_Y += [ED_data[n][J]]

        figData_QT_X = np.asarray(figData_QT_X); figData_QT_Y = np.asarray(figData_QT_Y); figData_QT_YErr = np.asarray(figData_QT_YErr)
        figData_ED_X = np.asarray(figData_ED_X); figData_ED_Y = np.asarray(figData_ED_Y)

        if LINEAR_FIT:
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

            X_fit = np.linspace(min(figData_QT_X), max(figData_QT_X), 10000)

            paramsED, cv = scipy.optimize.curve_fit(linFunc, figData_ED_X, figData_ED_Y, (1.0, 0.1))
            mED, bED = paramsED
            residuals = figData_ED_Y - linFunc(figData_ED_X, mED, bED)
            SSE = np.sum(residuals**2)
            diff = figData_ED_Y - figData_ED_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_ED = 1 - SSE/SST

            
            paramsQT, cv = scipy.optimize.curve_fit(linFunc, figData_QT_X, figData_QT_Y, (1.0, 0.1))
            mQT, bQT = paramsQT
            residuals = figData_QT_Y - linFunc(figData_QT_X, mQT, bQT)
            SSE = np.sum(residuals**2)
            diff = figData_QT_Y - figData_QT_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_QT = 1 - SSE/SST

            subfig1.plot(X_fit,  linFunc(X_fit, mED, bED), lw = 4, ls = "solid", markersize = 0, marker = "h", color = "blue", label = r"ED: $R^2 = %f$" % R2_ED)
            subfig1.plot(X_fit,  linFunc(X_fit, mQT, bQT), lw = 4, ls = "dashed", markersize = 0, marker = "s", color = "red", label = r"QT: $R^2 = %f$" % R2_QT)

            subfig1.plot(figData_ED_X, figData_ED_Y, lw = 0, ls = "solid", markersize = 16, marker = "h", color = "blue", label = "ED: Werte")
            if QT_WITH_ERROR: subfig1.errorbar(figData_QT_X, figData_QT_Y, yerr = figData_QT_YErr, lw = 0, ls = "dashed", markersize = 8, elinewidth = 3, fmt = "s", color = "red", label = "QT: Werte")
            else: subfig1.plot(figData_QT_X, figData_QT_Y, lw = 0, ls = "dashed", markersize = 8, marker = "s", color = "red", label = "QT: Werte")

            subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
            if J_SUM_SCALE: subfig1.set_ylabel(r'$\Delta_{SG}$ / ($J_1 + J_2$)', fontsize = labelfontsize)
            else: subfig1.set_ylabel(r'$\Delta_{SG}$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_title(r'Spingap Energien $\Delta_{SG}$ bei $J_1 / J_2$ = ' + J_val, fontsize = titlefontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

            fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_linear" + ".png")
            # fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_linear" + ".pdf")
            plt.close(fig1)

        if ONE_OVER_SQRTN_FIT:
            figData_QT_X_fit = np.asarray(figData_QT_X)
            figData_QT_X_fit = 1 / np.sqrt(figData_QT_X_fit)
            figData_ED_X_fit = np.asarray(figData_ED_X)
            figData_ED_X_fit = 1 / np.sqrt(figData_ED_X_fit)

            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

            #X_fit = np.linspace(min(figData_QT_X), max(figData_QT_X_fit), 10000)
            X_fit = np.linspace(0, max(figData_QT_X_fit), 10000)

            paramsED, cv = scipy.optimize.curve_fit(linFunc, figData_ED_X_fit, figData_ED_Y, (1.0, 0.1))
            mED, bED = paramsED
            residuals = figData_ED_Y - linFunc(figData_ED_X_fit, mED, bED)
            SSE = np.sum(residuals**2)
            diff = figData_ED_Y - figData_ED_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_ED = 1 - SSE/SST

            
            paramsQT, cv = scipy.optimize.curve_fit(linFunc, figData_QT_X_fit, figData_QT_Y, (1.0, 0.1))
            mQT, bQT = paramsQT
            residuals = figData_QT_Y - linFunc(figData_QT_X_fit, mQT, bQT)
            SSE = np.sum(residuals**2)
            diff = figData_QT_Y - figData_QT_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_QT = 1 - SSE/SST

            subfig1.plot(X_fit,  linFunc(X_fit, mED, bED), lw = 4, ls = "solid", markersize = 0, marker = "h", color = "blue", label = r"ED: $R^2 = %f$" % R2_ED)
            subfig1.plot(X_fit,  linFunc(X_fit, mQT, bQT), lw = 4, ls = "dashed", markersize = 0, marker = "s", color = "red", label = r"QT: $R^2 = %f$" % R2_QT)

            subfig1.plot(figData_ED_X_fit, figData_ED_Y, lw = 0, ls = "solid", markersize = 16, marker = "h", color = "blue", label = "ED")
            if QT_WITH_ERROR: subfig1.errorbar(figData_QT_X_fit, figData_QT_Y, yerr = figData_QT_YErr, lw = 0, ls = "dashed", markersize = 8, elinewidth = 3, fmt = "s", color = "red", label = "QT: Werte")
            else: subfig1.plot(figData_QT_X_fit, figData_QT_Y, lw = 0, ls = "dashed", markersize = 8, marker = "s", color = "red", label = "QT: Werte")

            subfig1.set_xlabel(r'$1/\sqrt{N}$', fontsize = labelfontsize)
            if J_SUM_SCALE: subfig1.set_ylabel(r'$\Delta_{SG}$ / ($J_1 + J_2$)', fontsize = labelfontsize)
            else: subfig1.set_ylabel(r'$\Delta_{SG}$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_title(r'Spingap Energien $\Delta_{SG}$ bei $J_1 / J_2$ = ' + J_val, fontsize = titlefontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

            fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_sqrt" + ".png")
            # fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_sqrt" + ".pdf")
            plt.close(fig1)

        if ONE_OVER_N_FIT:
            figData_QT_X_fit = np.asarray(figData_QT_X)
            figData_QT_X_fit = 1 / figData_QT_X_fit
            figData_ED_X_fit = np.asarray(figData_ED_X)
            figData_ED_X_fit = 1 / figData_ED_X_fit

            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

            #X_fit = np.linspace(min(figData_QT_X), max(figData_QT_X_fit), 10000)
            X_fit = np.linspace(0, max(figData_QT_X_fit), 10000)

            paramsED, cv = scipy.optimize.curve_fit(linFunc, figData_ED_X_fit, figData_ED_Y, (1.0, 0.1))
            mED, bED = paramsED
            residuals = figData_ED_Y - linFunc(figData_ED_X_fit, mED, bED)
            SSE = np.sum(residuals**2)
            diff = figData_ED_Y - figData_ED_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_ED = 1 - SSE/SST

            
            paramsQT, cv = scipy.optimize.curve_fit(linFunc, figData_QT_X_fit, figData_QT_Y, (1.0, 0.1))
            mQT, bQT = paramsQT
            residuals = figData_QT_Y - linFunc(figData_QT_X_fit, mQT, bQT)
            SSE = np.sum(residuals**2)
            diff = figData_QT_Y - figData_QT_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_QT = 1 - SSE/SST

            subfig1.plot(X_fit,  linFunc(X_fit, mED, bED), lw = 4, ls = "solid", markersize = 0, marker = "h", color = "blue", label = r"ED: $R^2 = %f$" % R2_ED)
            subfig1.plot(X_fit,  linFunc(X_fit, mQT, bQT), lw = 4, ls = "dashed", markersize = 0, marker = "s", color = "red", label = r"QT: $R^2 = %f$" % R2_QT)

            subfig1.plot(figData_ED_X_fit, figData_ED_Y, lw = 0, ls = "solid", markersize = 16, marker = "h", color = "blue", label = "ED")
            if QT_WITH_ERROR: subfig1.errorbar(figData_QT_X_fit, figData_QT_Y, yerr = figData_QT_YErr, lw = 0, ls = "dashed", markersize = 8, elinewidth = 3, fmt = "s", color = "red", label = "QT: Werte")
            else: subfig1.plot(figData_QT_X_fit, figData_QT_Y, lw = 0, ls = "dashed", markersize = 8, marker = "s", color = "red", label = "QT: Werte")

            subfig1.set_xlabel(r'$1/N$', fontsize = labelfontsize)
            if J_SUM_SCALE: subfig1.set_ylabel(r'$\Delta_{SG}$ / ($J_1 + J_2$)', fontsize = labelfontsize)
            else: subfig1.set_ylabel(r'$\Delta_{SG}$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_title(r'Spingap Energien $\Delta_{SG}$ bei $J_1 / J_2$ = ' + J_val, fontsize = titlefontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

            fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_N" + ".png")
            # fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_N" + ".pdf")
            plt.close(fig1)

        if ONE_OVER_N2_FIT:
            figData_QT_X_fit = np.asarray(figData_QT_X)
            figData_QT_X_fit = 1 / figData_QT_X_fit**2
            figData_ED_X_fit = np.asarray(figData_ED_X)
            figData_ED_X_fit = 1 / figData_ED_X_fit**2

            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

            #X_fit = np.linspace(min(figData_QT_X), max(figData_QT_X_fit), 10000)
            X_fit = np.linspace(0, max(figData_QT_X_fit), 10000)

            paramsED, cv = scipy.optimize.curve_fit(linFunc, figData_ED_X_fit, figData_ED_Y, (1.0, 0.1))
            mED, bED = paramsED
            residuals = figData_ED_Y - linFunc(figData_ED_X_fit, mED, bED)
            SSE = np.sum(residuals**2)
            diff = figData_ED_Y - figData_ED_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_ED = 1 - SSE/SST

            
            paramsQT, cv = scipy.optimize.curve_fit(linFunc, figData_QT_X_fit, figData_QT_Y, (1.0, 0.1))
            mQT, bQT = paramsQT
            residuals = figData_QT_Y - linFunc(figData_QT_X_fit, mQT, bQT)
            SSE = np.sum(residuals**2)
            diff = figData_QT_Y - figData_QT_Y.mean()
            square_diff = diff ** 2
            SST = square_diff.sum()
            R2_QT = 1 - SSE/SST

            subfig1.plot(X_fit,  linFunc(X_fit, mED, bED), lw = 4, ls = "solid", markersize = 0, marker = "h", color = "blue", label = r"ED: $R^2 = %f$" % R2_ED)
            subfig1.plot(X_fit,  linFunc(X_fit, mQT, bQT), lw = 4, ls = "dashed", markersize = 0, marker = "s", color = "red", label = r"QT: $R^2 = %f$" % R2_QT)

            subfig1.plot(figData_ED_X_fit, figData_ED_Y, lw = 0, ls = "solid", markersize = 16, marker = "h", color = "blue", label = "ED")
            if QT_WITH_ERROR: subfig1.errorbar(figData_QT_X_fit, figData_QT_Y, yerr = figData_QT_YErr, lw = 0, ls = "dashed", markersize = 8, elinewidth = 3, fmt = "s", color = "red", label = "QT: Werte")
            else: subfig1.plot(figData_QT_X_fit, figData_QT_Y, lw = 0, ls = "dashed", markersize = 8, marker = "s", color = "red", label = "QT: Werte")

            subfig1.set_xlabel(r'$1/N^2$', fontsize = labelfontsize)
            if J_SUM_SCALE: subfig1.set_ylabel(r'$\Delta_{SG}$ / ($J_1 + J_2$)', fontsize = labelfontsize)
            else: subfig1.set_ylabel(r'$\Delta_{SG}$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_title(r'Spingap Energien $\Delta_{SG}$ bei $J_1 / J_2$ = ' + J_val, fontsize = titlefontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)

            fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_N2" + ".png")
            # fig1.savefig("results/" + "SGExtrapFit/J_" + J_val + "_N2" + ".pdf")
            plt.close(fig1)

    end_time = time.time()
    print("done; this took %s" % format_time(start_time, end_time))
    print()
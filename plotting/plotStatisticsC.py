import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import sys
import os
plt.rcParams['text.usetex'] = True

N_color = []
N_color_LOW = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color_HIGH = [("18", "tomato"), ("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta"), ("28", "brown"), ("30", "purple"), ("32", "tomato")]
n_color = [("1", "red"), ("2", "blue"), ("3", "green"), ("4", "tomato")]
colors = ["red", "blue", "green", "magenta", "tomato", "brown", "purple"]

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

line_width = 4
marker_size = 0
alph = 0.1

valid_step_sizes = [0.1, 0.01, 0.001]
start = 0.0
ends = [0.1, 0.15, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0]

def get_set_size_index(stepsize):
    for i in range(len(valid_step_sizes)):
        if stepsize == valid_step_sizes[i]:
            return i
    else: return -1

def sort_step_sizes(X):
    length = len(X)
    for i in range(length-1):
        for j in range(0, length-i-1):
            x1 = X[j][0]; x2 = X[j+1][0]
            if x1 > x2:
                X[j], X[j+1] = X[j+1], X[j]
    return X

def plot_n_for_each_N():
    print("plotting dependance of n for fixed N (as error bands) ...")
    for N, NC in N_color:
        print("N = " + N)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        # ED results
        file = open("results/" + N + "/data/data_specific_heat_J_const.txt", 'r')
        lines = file.readlines()
        linesJ = lines[0][len("J1/J2: "):-1]
        #linesh = lines[2][len("h: "):-1]
        lbl = "ED: N = " + N
        X = []
        Y = []
        for i in range(9,len(lines)):
            x, y = lines[i].split("\t")
            if float(x) <= 0: continue
            X += [float(x)]
            Y += [float(y)]
        subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "black", label = lbl)
        # QT results
        for n, nc in n_color:
            file = open("results/" + N + "/data/" + n +"_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            linesh = lines[2][len("h: "):-1]
            lbl = "n = " + n
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                if float(x) <= 0: continue
                X += [float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            subfig1.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = nc, label = "QT: n = " + n, alpha = 1.0)
            X = np.asarray(X)
            Y = np.asarray(Y)
            YErr = np.asarray(YErr)
            subfig1.fill_between(X, Y - YErr, Y + YErr, color = nc, alpha = 0.1)
        # saving
        subfig1.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
        subfig1.set_ylabel(r'$C/N$ in $J_2$', fontsize = 25)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        subfig1.set_title(r"$C/N$ für N = " + N + r" mit $J_1$ / $J_2$ = " + linesJ, fontsize = 25)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)
        filename = "N_" + N + "_n"
        for n, nc in n_color:
            filename += "_" + n
        for end in ends:
            subfig1.set_xlim(start, end)
            plt.savefig("results/QT_stats/C_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.savefig("results/QT_stats/C_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig1)

def plot_n_for_each_N_sigma_abs(): 
    print("plotting dependance of n for fixed N (plotting only sigma abs) ...")
    for N, NC in N_color:
        print("N = " + N)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        x_min = 42069
        # QT results
        for n, nc in n_color:
            file = open("results/" + N + "/data/" + n +"_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            #linesh = lines[2][len("h: "):-1]
            #lbl = "n = " + n
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                if float(x) <= 0: continue
                X += [1.0/float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = nc, label = "QT: n = " + n)
            if min(X) < x_min: x_min = min(X)
        # saving
        subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\sigma_{abs}$ / $J_2$', fontsize = labelfontsize)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        subfig1.set_title(r"$\sigma_{abs}$ von $C/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + linesJ, fontsize = titlefontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        filename = "N_" + N + "_n"
        for n, nc in n_color:
            filename += "_" + n
        for end in ends:
            subfig1.set_xlim(x_min, end)
            plt.savefig("results/QT_stats/C_sigma_abs_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.savefig("results/QT_stats/C_sigma_abs_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig1)

def plot_n_for_each_N_sigma_rel(): 
    print("plotting dependance of n for fixed N (plotting only sigma rel) ...")
    for N, NC in N_color:
        print("N = " + N)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        x_min = 42069
        # QT results
        for n, nc in n_color:
            file = open("results/" + N + "/data/" + n +"_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            #linesh = lines[2][len("h: "):-1]
            #lbl = "n = " + n
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                if float(x) <= 0: continue
                if not float(y) > 0.0: continue
                X += [1.0/float(x)]
                Y += [float(y)]
                YErr += [float(yErr)/float(y)]
            subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = nc, label = "QT: n = " + n)
            if min(X) < x_min: x_min = min(X)
        # saving
        subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\sigma_{rel}$ / $J_2$', fontsize = labelfontsize)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        subfig1.set_title(r"$\sigma_{rel}$ von $C/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + linesJ, fontsize = titlefontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        filename = "N_" + N + "_n"
        for n, nc in n_color:
            filename += "_" + n
        for end in ends:
            subfig1.set_xlim(x_min, end)
            plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig1)

def plot_N_for_each_n():
    print("plotting dependance of N for fixed n (as error bands) ...")
    for n, nc in n_color:
        print("n = " + n)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        for N, NC in N_color:
            # ED results
            file = open("results/" + N + "/data/data_specific_heat_J_const.txt", 'r')
            lines = file.readlines()
            linesJ = lines[0][len("J1/J2: "):-1]
            #linesh = lines[2][len("h: "):-1]
            lbl = " ED: N = " + N
            X = []
            Y = []
            for i in range(9,len(lines)):
                x, y = lines[i].split("\t")
                if float(x) <= 0: continue
                X += [float(x)]
                Y += [float(y)]
            subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 1, marker = "o", color = NC, label = lbl)
            # QT results
            file = open("results/" + N + "/data/" + n +"_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            #linesh = lines[2][len("h: "):-1]
            lbl = "n = " + n
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                if float(x) <= 0: continue
                X += [float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            subfig1.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = NC, label = "QT: N = " + N, alpha = 1.0)
            X = np.asarray(X)
            Y = np.asarray(Y)
            YErr = np.asarray(YErr)
            subfig1.fill_between(X, Y - YErr, Y + YErr, color = NC, alpha = 0.1)
        # saving
        subfig1.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
        subfig1.set_ylabel(r'$C/N$ in $J_2$', fontsize = 25)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        subfig1.set_title(r"$C/N$ bei der QT mit $J_1$ / $J_2$ = " + linesJ + " und " + n + " Startvektoren", fontsize = 25)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)
        filename = "n_" + n + "_N"
        for N, NC in N_color:
            filename += "_" + N
        for end in ends:
            subfig1.set_xlim(start, end)
            plt.savefig("results/QT_stats/C_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.savefig("results/QT_stats/C_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig1)

def plot_N_for_each_n_sigma_abs():
    print("plotting dependance of N for fixed n (plotting only sigma abs) ...")
    for n, nc in n_color:
        print("n = " + n)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        x_min = 42069
        # QT results
        for N, NC in N_color:
            # QT results
            file = open("results/" + N + "/data/" + n +"_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            #linesh = lines[2][len("h: "):-1]
            #lbl = "n = " + n
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                if float(x) <= 0: continue
                X += [1.0/float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = NC, label = "QT: N = " + N, alpha = 1.0)
            if min(X) < x_min: x_min = min(X)
        # saving
        subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\sigma_{abs}$ / $J_2$', fontsize = labelfontsize)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        vec = "einem Startvektor"
        if int(n) > 1: vec = n + " Startvektoren"
        subfig1.set_title(r"$\sigma_{abs}$ von $C/N$ mit $J_1$ / $J_2$ = " + linesJ + " und " + vec, fontsize = titlefontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        filename = "n_" + n + "_N"
        for N, NC in N_color:
            filename += "_" + N
        for end in ends:
            subfig1.set_xlim(x_min, end)
            plt.savefig("results/QT_stats/C_sigma_abs_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.savefig("results/QT_stats/C_sigma_abs_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig1)

def plot_N_for_each_n_sigma_rel():
    print("plotting dependance of N for fixed n (plotting only sigma rel) ...")
    for n, nc in n_color:
        print("n = " + n)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        x_min = 42069
        # QT results
        for N, NC in N_color:
            # QT results
            file = open("results/" + N + "/data/" + n +"_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            #linesh = lines[2][len("h: "):-1]
            #lbl = "n = " + n
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                if float(x) <= 0: continue
                if not float(y) > 0.0: continue
                X += [1.0/float(x)]
                Y += [float(y)]
                YErr += [float(yErr)/float(y)]
            subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = NC, label = "QT: N = " + N, alpha = 1.0)
            if min(X) < x_min: x_min = min(X)
        # saving
        subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$\sigma_{rel}$ / $J_2$', fontsize = labelfontsize)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        vec = "einem Startvektor"
        if int(n) > 1: vec = n + " Startvektoren"
        subfig1.set_title(r"$\sigma_{rel}$ von $C/N$ mit $J_1$ / $J_2$ = " + linesJ + " und " + vec, fontsize = titlefontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        filename = "n_" + n + "_N"
        for N, NC in N_color:
            filename += "_" + N
        for end in ends:
            subfig1.set_xlim(x_min, end)
            plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig1)

def plot_delta_ED():
    print("plotting \Delta E ED - QT for fixed N ...")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    used_N = "N"
    x_min = 0
    for N, NC in N_color:
        print("N = " + N)
        # QT results
        fileQT = open("results/" + N + "/data/1_data_specific_heat_J_const_QT.txt", 'r')
        linesQT = fileQT.readlines()
        linesJQT = linesQT[1][len("J1/J2: "):-1]
        #linesh = lines[2][len("h: "):-1]
        #lbl = "n = " + n
        X_QT = []
        Y_QT = []
        for i in range(7,len(linesQT)):
            x, y, yErr = linesQT[i].split("\t")
            if float(x) <= 0: continue
            X_QT += [1.0/float(x)]
            Y_QT += [float(y)]
        # ED results
        fileED = open("results/" + N + "/data/data_specific_heat_J_const.txt", 'r')
        linesED = fileED.readlines()
        linesJED = linesED[1][len("J1/J2: "):-1]
        #linesh = lines[2][len("h: "):-1]
        #lbl = "n = " + n
        X_ED = []
        Y_ED = []
        for i in range(9,len(linesED)):
            x, y = linesED[i].split("\t")
            if float(x) <= 0: continue
            X_ED += [1.0/float(x)]
            Y_ED += [float(y)]

        X = []
        Y_abs = []; Y_rel = []

        # print(len(X_QT))
        # print(len(Y_QT))
        # print(len(X_ED))
        # print(len(Y_ED))

        for i in range(0,len(X_QT)):
            try:
                if X_QT[i] > 0.5: continue
                if X_QT[i] != X_ED[i]: continue
                Y_abs += [abs(Y_QT[i] - Y_ED[i])]
                Y_rel += [abs(Y_QT[i] - Y_ED[i]) / Y_ED[i]]
                X += [X_QT[i]]
            except:
                pass

        if min(X) > x_min: x_min = min(X)

        subfig1.plot(X, Y_abs, lw = line_width, ls = "solid", markersize = 0, marker = "o", color = NC, label = r"N = " + N, alpha = 1.0) # label = r"$\vert$QT-ED$\vert$: N = " + N
        subfig1.plot(X, Y_rel, lw = line_width, ls = "dashed", markersize = 0, marker = "o", color = NC)

        used_N += "_" + N

    # saving
    subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'$\Delta_{QT-ED}$ / $J_2$', fontsize = labelfontsize)
    #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
    subfig1.set_title(r"$\Delta_{QT-ED}$ von $C/N$ bei $J_1$ / $J_2$ = " + linesJQT, fontsize = titlefontsize)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    for end in ends:
        subfig1.set_xlim(x_min, end)
        plt.savefig("results/QT_stats/C_delta_ED_QT_" + used_N + "_J" + linesJQT + "_" + str(start) + "_" + str(end) + ".png")
        plt.savefig("results/QT_stats/C_delta_ED_QT_" + used_N + "_J" + linesJQT + "_" + str(start) + "_" + str(end) + ".pdf")
    plt.close(fig1)

def plot_step_size():
    print("plotting step size dependance for fixed N ...")
    for N, NC in N_color:
        print("N = " + N)
        for n, nc in n_color:
            if n != "1": continue
            print("n = " + n)
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            # ED results
            file = open("results/" + N + "/data/data_specific_heat_J_const.txt", 'r')
            lines = file.readlines()
            linesJ = lines[0][len("J1/J2: "):-1]
            lbl = "ED: N = " + N
            X = []
            Y = []
            for i in range(9,len(lines)):
                x, y = lines[i].split("\t")
                if float(x) <= 0: continue
                X += [1.0/float(x)]
                Y += [float(y)]
            file.close()
            subfig1.plot(X, Y, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = "black")#, label = lbl)
            # QT results
            filenum = 0
            used_step_sizes = ""
            step_sizes = []
            for filename in os.listdir("results/" + N + "/data/step_size_data/"):
                if n + "_data_specific_heat_J_const_QT_step" in filename:# and "_data_specific_heat_J_const_QT_step" not in filename:
                    stepsize = filename[len(n + "_data_specific_heat_J_const_QT_step"): - len(".txt")]
                    if float(stepsize) not in valid_step_sizes: continue
                    filenum = get_set_size_index(float(stepsize))
                    step_sizes += [(float(stepsize), colors[filenum])]
                    file = open("results/" + N + "/data/step_size_data/" + n + "_data_specific_heat_J_const_QT_step" + stepsize + ".txt", 'r')
                    lines = file.readlines()
                    X = []; Y = []; YErr = []
                    for i in range(7,len(lines)):
                        x, y, yErr = lines[i].split("\t")
                        if float(x) <= 0: continue
                        X += [1.0/float(x)]
                        Y += [float(y)]
                        YErr += [float(yErr)]
                    file.close()
                    subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[filenum])#, label = "QT: " + str(float(stepsize)))
                    X = np.asarray(X)
                    Y = np.asarray(Y)
                    YErr = np.asarray(YErr)
                    subfig1.fill_between(X, Y - YErr, Y + YErr, color = colors[filenum], alpha = alph)
                    filenum += 1
                    used_step_sizes += "_" + str(float(stepsize))
            # saving
            subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_ylabel(r'$C/N$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_title(r"spezifische Wärmekapazität pro Spin $C/N$" + "bei unterschiedlichen\nSchrittweiten und 12 Mittelungen über einen Startvektor", fontsize = titlefontsize)
            subfig1.axhline(0, color = "grey")
            legend = []
            color_count = 0
            legend += [Line2D([0], [0], label = "ED: N = " + N, color = "black", ls = "solid", lw = line_width)]
            step_sizes = sort_step_sizes(step_sizes)
            for sz, c in step_sizes:
                legend += [Line2D([0], [0], label = "QT: " + str(sz), color = c, ls = "dashed", lw = line_width)]
            handles, labels = plt.gca().get_legend_handles_labels()
            handles.extend(legend)
            subfig1.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)
            subfig1.tick_params(axis = "both", which = "major", labelsize = axisfontsize)
            for end in ends:
                subfig1.set_xlim(start, end)
                plt.savefig("results/QT_stats/C_N_" + N + "_n_" + n + "_step_size" + used_step_sizes + "_J_" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/C_N_" + N + "_n_" + n + "_step_size" + used_step_sizes + "_J_" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.close(fig1)
        print()

if __name__ == "__main__":
    print("plotting specific heat (constant J1/J2, funtion of T):")
    N_color = [("10", "red"), ("12", "blue"), ("14", "green"), ("16", "purple"), ("18", "tomato")] # ("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), 
    # start = float(sys.argv[1])
    # end = float(sys.argv[2])
    # regime = sys.argv[3]

    # if regime == "low": N_color = N_color_LOW
    # elif regime == "high": N_color = N_color_HIGH
    # else: exit()

    # plot_n_for_each_N()
    # print()
    # plot_n_for_each_N_sigma_abs()
    # print()
    # plot_n_for_each_N_sigma_rel()
    # print()
    # plot_N_for_each_n()
    # print()
    # N_color = [("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
    # plot_N_for_each_n_sigma_abs()
    # print()
    # plot_N_for_each_n_sigma_rel()
    # print()
    # N_color_HIGH = [("18", "red"), ("20", "blue"), ("22", "green")]
    plot_delta_ED()
    print()
    # plot_step_size()
    # print()

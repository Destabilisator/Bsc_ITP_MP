import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import sys
import os
import gc
plt.rcParams['text.usetex'] = True

N_color = []
N_color_LOW = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color_HIGH = [("18", "tomato"), ("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta"), ("28", "brown"), ("30", "purple"), ("32", "tomato")]
n_color = [("1", "red"), ("2", "blue"), ("3", "green"), ("4", "tomato")]
colors = ["red", "blue", "green", "magenta", "tomato", "brown", "purple", "cyan"]
N_Array = []

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

line_width = 4
marker_size = 0
alph = 0.1
max_n = 5

valid_step_sizes = [0.1, 0.01, 0.001]
start = 0.0
ends = [0.1, 0.15, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0]
valid_J = [0.010000, 0.358600, 1.006000, 1.753000]

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
        if abs(J-j) <= 0.0001:
            return True
    return False

def Npos(N):
    for i in range(len(N_Array)):
        if N_Array[i] == N:
            return i
    return -1

def C_plot_n_for_each_N_sigma_rel(): 
    print("C: plotting dependance of n for fixed N (plotting only sigma rel) ...")
    s = -1
    # try:
    for N in N_Array:
        N = str(N)
        # QT results
        J_cnt = 0
        for filename in os.listdir("results/" + N + "/data/excitation_energies_data/1/"):
            x_min = 42069
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
                    if not float(y) > 0.0: continue
                    X += [1.0/float(x)]; X_arr[n] += [1.0/float(x)]
                    Y += [float(y)]; Y_arr[n] += [float(y)]
                file.close()
                if min(X) < x_min: x_min = min(X)
            
            X = X_arr[0]
            filename = "N_" + N + "_n"
            color_count = 0
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            for samp in range(1, int(max_n / 2) + 1):
                # for n in range(max_n):
                #     subfig1.plot(X, Y_arr[n], lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = "black")
                s = samp
                print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                Y, YErr = avgData(samp, Y_arr)
                Y = np.asarray(Y); YErr = np.asarray(YErr)
                YErr = YErr / Y
                # ploting
                subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "n = %s" % str(samp))
                color_count += 1
                # subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count], label = r"$n = %s$" % str(samp))
                # color_count += 1
                if color_count >= len(colors): break
                # saving
                subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
                subfig1.set_ylabel(r'$\sigma_{rel}$ / $J_2$', fontsize = labelfontsize)
                subfig1.set_title(r"$\sigma_{rel}$ von $C/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
                subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
                filename += "_" + str(samp)
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...               ")

def C_plot_N_for_each_n_sigma_rel():
    print("C: plotting dependance of N for fixed n (plotting only sigma rel) ...")
    s = -1
    # try:
    for filename in os.listdir("results/" + str(N_Array[0]) + "/data/excitation_energies_data/1/"):
        x_min = 42069
        if ".png" in filename or ".pdf" in filename: continue
        if "ED" in filename: continue
        if "placeholder" in filename: continue
        J = filename[len("C_J"):-len("QT.txt")]
        if not isValid(float(J)): continue
        N_SAMP = [[]] * len(N_Array)
        for N in N_Array:
            # QT results
            J_cnt = 0
            for filename in os.listdir("results/" + str(N) + "/data/excitation_energies_data/1/"):
                x_min = 42069
                if ".png" in filename or ".pdf" in filename: continue
                if "ED" in filename: continue
                if "placeholder" in filename: continue
                J = filename[len("C_J"):-len("QT.txt")]
                if not isValid(float(J)): continue
                J_cnt += 1
                print("get data: N: %s, J: %s, %i/%i          \r" % (N, str(J), J_cnt, len(valid_J)), end = "")
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
                        if not float(y) > 0.0: continue
                        X += [1.0/float(x)]; X_arr[n] += [1.0/float(x)]
                        Y += [float(y)]; Y_arr[n] += [float(y)]
                    file.close()
                    if min(X) < x_min: x_min = min(X)
                
                X = X_arr[0]
                color_count = 0
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                sampArr = [[]]
                for samp in range(1, int(max_n / 2) + 1):
                    print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                    Y, YErr = avgData(samp, Y_arr)
                    Y = np.asarray(Y); YErr = np.asarray(YErr)
                    YErr = YErr / Y
                    sampArr += [YErr]
                N_SAMP += [sampArr]
        for samp in range(len(N_SAMP[0])):
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            filename = "n_" + str(samp) + "_N"
            print("plotting: %i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
            for N in range(len(N_Array)):
                YErr = N_SAMP[N][samp]
                # ploting
                subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "N = " + str(N_Array[0]))
                color_count += 1
                filename += "_" + str(N_Array[0])
                # subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count], label = r"$n = %s$" % str(samp))
                # color_count += 1
                if color_count >= len(colors): break
            # saving
            subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_ylabel(r'$\sigma_{rel}$ / $J_2$', fontsize = labelfontsize)
            vec = "einem Startvektor"
            if int(samp) > 1: vec = str(samp) + " Startvektoren"
            subfig1.set_title(r"$\sigma_{rel}$ von $C/N$ bei der QT mit Mittelung über" + vec + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
            subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/C_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...               ")

def X_plot_n_for_each_N_sigma_rel(): 
    print("X: plotting dependance of n for fixed N (plotting only sigma rel) ...")
    s = -1
    try: 
        for N in N_Array:
            N = str(N)
            # QT results
            J_cnt = 0
            for filename in os.listdir("results/" + N + "/data/spin_gap_data/1/"):
                x_min = 42069
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
                        if not float(y) > 0.0: continue
                        X += [1.0/float(x)]; X_arr[n] += [1/float(x)]
                        Y += [float(y)]; Y_arr[n] += [float(y)]
                    file.close()
                    if min(X) < x_min: x_min = min(X)
                
                X = X_arr[0]
                filename = "N_" + N + "_n"
                color_count = 0
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                for samp in range(1, int(max_n / 2) + 1):
                    # for n in range(max_n):
                    # subfig1.plot(X, Y_arr[n], lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = "black")
                    s = samp
                    print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                    Y, YErr = avgData(samp, Y_arr)
                    Y = np.asarray(Y); YErr = np.asarray(YErr)
                    YErr = YErr / Y
                    # ploting
                    subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "n = %s" % str(samp))
                    color_count += 1
                    # subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count], label = r"$n = %s$" % str(samp))
                    # color_count += 1
                    if color_count >= len(colors): break
                    # saving
                    subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
                    subfig1.set_ylabel(r'$\sigma_{rel}$ / $J_2$', fontsize = labelfontsize)
                    subfig1.set_title(r"$\sigma_{rel}$" + " von $\\chi/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
                    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
                    subfig1.axhline(0, color = "grey")
                    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
                    filename += "_" + str(samp)
                for end in ends:
                    subfig1.set_xlim(x_min, end)
                    plt.savefig("results/QT_stats/X_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                    plt.savefig("results/QT_stats/X_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
                plt.cla()
                plt.clf()
                plt.close(fig1)
                plt.close('all')
                del fig1, subfig1
                gc.collect()
    except:
        print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...               ")

def X_plot_N_for_each_n_sigma_rel():
    print("X: plotting dependance of N for fixed n (plotting only sigma rel) ...")
    s = -1
    # try:
    for filename in os.listdir("results/" + str(N_Array[0]) + "/data/spin_gap_data/1/"):
        x_min = 42069
        if ".png" in filename or ".pdf" in filename: continue
        if "ED" in filename: continue
        if "placeholder" in filename: continue
        J = filename[len("C_J"):-len("QT.txt")]
        if not isValid(float(J)): continue
        N_SAMP = [[]] * len(N_Array)
        for N in N_Array:
            # QT results
            J_cnt = 0
            for filename in os.listdir("results/" + str(N) + "/data/spin_gap_data/1/"):
                x_min = 42069
                if ".png" in filename or ".pdf" in filename: continue
                if "ED" in filename: continue
                if "placeholder" in filename: continue
                J = filename[len("C_J"):-len("QT.txt")]
                if not isValid(float(J)): continue
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
                        if not float(y) > 0.0: continue
                        X += [1.0/float(x)]; X_arr[n] += [1.0/float(x)]
                        Y += [float(y)]; Y_arr[n] += [float(y)]
                    file.close()
                    if min(X) < x_min: x_min = min(X)
                
                X = X_arr[0]
                color_count = 0
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                sampArr = [[]]
                for samp in range(1, int(max_n / 2) + 1):
                    print("%i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
                    Y, YErr = avgData(samp, Y_arr)
                    Y = np.asarray(Y); YErr = np.asarray(YErr)
                    YErr = YErr / Y
                    sampArr += [YErr]
                N_SAMP += [sampArr]
        for samp in range(len(N_SAMP[0])):
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            filename = "n_" + str(samp) + "_N"
            print("plotting: %i sample: N: %s, J: %s, %i/%i          \r" % (samp, N, str(J), J_cnt, len(valid_J)), end = "")
            for N in range(len(N_Array)):
                YErr = N_SAMP[N][samp]
                # ploting
                subfig1.plot(X, YErr, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "N = " + str(N_Array[0]))
                color_count += 1
                filename += "_" + str(N_Array[0])
                # subfig1.plot(X, Y, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count], label = r"$n = %s$" % str(samp))
                # color_count += 1
                if color_count >= len(colors): break
            # saving
            subfig1.set_xlabel(r'$k_B T$ / $J_2$', fontsize = labelfontsize)
            subfig1.set_ylabel(r'$\sigma_{rel}$ / $J_2$', fontsize = labelfontsize)
            vec = "einem Startvektor"
            if int(samp) > 1: vec = str(samp) + " Startvektoren"
            subfig1.set_title(r"$\sigma_{rel}$" + " von $\\chi/N$ bei der QT mit Mittelung über" + vec + r" und $J_1$ / $J_2$ = " + str(J), fontsize = titlefontsize)
            subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            for end in ends:
                subfig1.set_xlim(x_min, end)
                plt.savefig("results/QT_stats/X_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".png")
                plt.savefig("results/QT_stats/X_sigma_rel_" + filename + "_J" + str(J) + "_" + str(start) + "_" + str(end) + ".pdf")
            plt.cla()
            plt.clf()
            plt.close(fig1)
            plt.close('all')
            del fig1, subfig1
            gc.collect()
    # except:
    #     print("%i sample: N: %s, J: %s, %i/%i: FAILED" % (s, N, str(J), J_cnt, len(valid_J)))
    print("done...               ")

if __name__ == "__main__":
    N_Array = [16]
    C_plot_n_for_each_N_sigma_rel()
    print()
    C_plot_N_for_each_n_sigma_rel()
    print()
    X_plot_n_for_each_N_sigma_rel()
    print()
    X_plot_N_for_each_n_sigma_rel()
    print()

import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('TkAgg')
import os
import sys
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy.optimize
# from typing import Tuple
import gc
plt.rcParams['text.usetex'] = True

colors = ["red", "blue", "green", "tomato", "purple"]

line_width = 3
marker_size = 5

counter = 0

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

timestepfontsize = 20

SAVE_FULL_PLOTS = True
ONLY_NEWEST_RUN = False

samples = ["1"]#, "2", "3"]
stepsizes = ["0.100000"] # ["0.100000", "0.010000"]
cores = ["1"] # ["1", "2", "5", "10"]

def strip_data(N,T):
    arr = [0] * 35
    out_N = []; out_T = []
    for i in range(len(N)):
        arr[N[i]] = T[i]
    for n in range(len(arr)):
        if arr[n] != 0:
            out_N += [n]; out_T += [arr[n]]
    return out_N, out_T
    
def sort_data(N, T):
    if ONLY_NEWEST_RUN: N, T = strip_data(N, T)
    length = len(N)
    for i in range(length-1):
        for j in range(0, length-i-1):
            N1 = N[j]; N2 = N[j+1]
            if N1 > N2:
                N[j], N[j+1] = N[j+1], N[j]
                T[j], T[j+1] = T[j+1], T[j]
    return N, T

def extrap_func(x: float, A: float, k: float, x_0: float, b: float) -> float:
    return A * np.exp(k * x - x_0) + b

def extrapolate_data(N, T, title):
    Y_fit = np.log(T)
    N = np.asarray(N)
    params = poly.polyfit(N, Y_fit, 1, w = N)# + 1/2 * N**2)#  + 1/6 * N**3) # , w = N + 1/2 * N**2
    fit_func = poly.Polynomial(params)
    X = np.linspace(6, 32, 1000)
    Y = fit_func(X)
    Y = np.exp(Y)

    if SAVE_FULL_PLOTS:
        global counter
        fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
        subfig2.plot(N, T, lw = 0, ls = "solid", markersize = 5, marker = "o", color = "black")
        subfig2.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = "black")
        subfig2.set_title(title, fontsize = 40)
        subfig2.set_yscale('log')
        fig2.savefig("./results/benchmarking/temp/" + str(counter) + ".png")
        counter += 1
        plt.close(fig2)

    return X, Y

def add_time_steps(subfig):
    # Namen aus DIN 1355-1
    width = 2; c = "black"; alph = 0.5
    # 1 second
    subfig.hlines(y = 1, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7, "1 Sek.", fontsize = timestepfontsize)
    # 1 minute
    subfig.hlines(y = 60, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7 * 60, "1 Min.", fontsize = timestepfontsize)
    # 1 hour
    subfig.hlines(y = 60 * 60, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7 * 60 * 60, "1 Std.", fontsize = timestepfontsize)
    # 1 day
    subfig.hlines(y = 60 * 60 * 24, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7 * 60 * 60 * 24, "1 Tg.", fontsize = timestepfontsize)
    # 1 week
    subfig.hlines(y = 60 * 60 * 24 * 7, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7 * 60 * 60 * 24 * 7, "1 Wo.", fontsize = timestepfontsize)
    # 1 month
    subfig.hlines(y = 60 * 60 * 24 * 31, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7 * 60 * 60 * 24 * 31, "31 Tg.", fontsize = timestepfontsize) # 1 Mon.
    # 1 year
    subfig.hlines(y = 60 * 60 * 24 * 365, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(4.8, 0.7 * 60 * 60 * 24 * 365, "1 J.", fontsize = timestepfontsize)

def add_size_steps(subfig):
    width = 2; c = "black"; alph = 0.5
    # 1 KB
    # subfig.hlines(y = 1024 ** 1 / 16, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    # subfig.text(5.2, 0.7 * 1024 ** 1 / 16, "1 KB", fontsize = timestepfontsize)
    # 1 MB
    subfig.hlines(y = 1024 ** 2 / 16, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(5.2, 0.7 * 1024 ** 2 / 16, "1 MB", fontsize = timestepfontsize)
    # 1 GB
    subfig.hlines(y = 1024 ** 3 / 16, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(5.2, 0.7 * 1024 ** 3 / 16, "1 GB", fontsize = timestepfontsize)
    # 1 TB
    subfig.hlines(y = 1024 ** 4 / 16, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(5.2, 0.7 * 1024 ** 4 / 16, "1 TB", fontsize = timestepfontsize)
    # 1 PB
    subfig.hlines(y = 1024 ** 5 / 16, xmin = 6.5, xmax = 32, lw = width, color = c, ls = "dashed", alpha = alph)
    subfig.text(5.2, 0.7 * 1024 ** 5 / 16, "1 PB", fontsize = timestepfontsize)

##### run time, full matrix #####
def RT_plot_raw_files():
    print("plotting raw files")
    for core in cores:
        # ED exac data
        ED_SG_file = open("./results/benchmarking/runtime/data/" + "ED_SG_cores_" + core + ".txt")
        ED_SG_lines = ED_SG_file.readlines()
        ED_SG_file.close()
        N_ED_SG = []; T_ED_SG = []
        for line in ED_SG_lines:
            n, t = line.split("\t")
            N_ED_SG += [int(n)]; T_ED_SG += [float(t)]
        N_ED_SG, T_ED_SG = sort_data(N_ED_SG, T_ED_SG)
        N_ED_SG_extrap, T_ED_SG_extrap = extrapolate_data(N_ED_SG, T_ED_SG, "ed sg full")

        for stepsize in stepsizes:
            # ED fit data
            ED_MJ_file = open("./results/benchmarking/runtime/data/" + "ED_MJ_step_" + stepsize + "_cores_" + core + ".txt")
            ED_MJ_lines = ED_MJ_file.readlines()
            ED_MJ_file.close()
            N_ED_MJ = []; T_ED_MJ = []
            for line in ED_MJ_lines:
                n, t = line.split("\t")
                N_ED_MJ += [int(n)]; T_ED_MJ += [float(t)]
            N_ED_MJ, T_ED_MJ = sort_data(N_ED_MJ, T_ED_MJ)
            N_ED_MJ_extrap, T_ED_MJ_extrap = extrapolate_data(N_ED_MJ, T_ED_MJ, "ed mj")

            for sample in samples:
                # QT data
                QT_file = open("./results/benchmarking/runtime/data/" + "QT_SG_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".txt")
                QT_lines = QT_file.readlines()
                QT_file.close()
                N_QT = []; T_QT = []
                for line in QT_lines:
                    n, t = line.split("\t")
                    N_QT += [int(n)]; T_QT += [float(t)]
                N_QT, T_QT = sort_data(N_QT, T_QT)
                N_QT_extrap, T_QT_extrap = extrapolate_data(N_QT, T_QT, "qt full")

                # plotting
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                color_count = 0
                subfig1.plot(N_ED_SG, T_ED_SG, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_ED_SG_extrap, T_ED_SG_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED-ev")
                color_count += 1
                subfig1.plot(N_ED_MJ, T_ED_MJ, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_ED_MJ_extrap, T_ED_MJ_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED-fit")
                color_count += 1
                subfig1.plot(N_QT, T_QT, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_QT_extrap, T_QT_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "QT-fit")
                color_count += 1
                subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
                subfig1.set_ylabel(r'$t$ in $s$', fontsize = labelfontsize)
                if sample == "1": vec = "Startvektor"
                else: vec = "Startvektoren"
                subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$" + "\n" + "Schrittweite (nur Fits) = %.2f, %s %s (nur QT)" % (float(stepsize), sample, vec), fontsize = titlefontsize)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
                # plt.xticks(fontsize = 25)
                # plt.yticks(fontsize = 25)
                subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
                add_time_steps(subfig1)
                fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_SG_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
                fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_SG_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".pdf")
                subfig1.set_yscale('log')
                fig1.savefig("./results/benchmarking/runtime/" + "log_QT_ED_SG_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
                fig1.savefig("./results/benchmarking/runtime/" + "log_QT_ED_SG_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".pdf")
                plt.close(fig1)

def RT_plot_only_ED():
    print("comparing only ED")
    for core in cores:
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        color_count = 0
        # ED exac data
        ED_SG_file = open("./results/benchmarking/runtime/data/" + "ED_SG_cores_" + core + ".txt")
        ED_SG_lines = ED_SG_file.readlines()
        ED_SG_file.close()
        N_ED_SG = []; T_ED_SG = []
        for line in ED_SG_lines:
            n, t = line.split("\t")
            N_ED_SG += [int(n)]; T_ED_SG += [float(t)]
        N_ED_SG, T_ED_SG = sort_data(N_ED_SG, T_ED_SG)
        N_ED_SG_extrap, T_ED_SG_extrap = extrapolate_data(N_ED_SG, T_ED_SG, "ed sg ed")

        subfig1.plot(N_ED_SG, T_ED_SG, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
        subfig1.plot(N_ED_SG_extrap, T_ED_SG_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED-ev")
        color_count += 1

        for stepsize in stepsizes:
            # ED fit data
            ED_MJ_file = open("./results/benchmarking/runtime/data/" + "ED_MJ_step_" + stepsize + "_cores_" + core + ".txt")
            ED_MJ_lines = ED_MJ_file.readlines()
            ED_MJ_file.close()
            N_ED_MJ = []; T_ED_MJ = []
            for line in ED_MJ_lines:
                n, t = line.split("\t")
                N_ED_MJ += [int(n)]; T_ED_MJ += [float(t)]
            N_ED_MJ, T_ED_MJ = sort_data(N_ED_MJ, T_ED_MJ)#, "ed mj")
            N_ED_MJ_extrap, T_ED_MJ_extrap = extrapolate_data(N_ED_MJ, T_ED_MJ, "ed mj ed")

            subfig1.plot(N_ED_MJ, T_ED_MJ, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
            subfig1.plot(N_ED_MJ_extrap, T_ED_MJ_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED-fit: Schrittweite %.2f" % float(stepsize))
            color_count += 1

        subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = labelfontsize)
        subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$", fontsize = titlefontsize)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        add_time_steps(subfig1)
        fig1.savefig("./results/benchmarking/runtime/" + "ED_SG_MJ" + "_cores_" + core + ".png")
        fig1.savefig("./results/benchmarking/runtime/" + "ED_SG_MJ" + "_cores_" + core + ".pdf")
        subfig1.set_yscale('log')
        fig1.savefig("./results/benchmarking/runtime/" + "log_ED_SG_MJ" + "_cores_" + core + ".png")
        fig1.savefig("./results/benchmarking/runtime/" + "log_ED_SG_MJ" + "_cores_" + core + ".pdf")
        plt.close(fig1)

def RT_plot_step_size_influence():
    print("comparing influence of stepsize")
    for core in cores:
        for sample in samples:
            # fig for stepsize influence
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            color_count = 0
            for stepsize in stepsizes:
                # ED fit data
                ED_MJ_file = open("./results/benchmarking/runtime/data/" + "ED_MJ_step_" + stepsize + "_cores_" + core + ".txt")
                ED_MJ_lines = ED_MJ_file.readlines()
                ED_MJ_file.close()
                N_ED_MJ = []; T_ED_MJ = []
                for line in ED_MJ_lines:
                    n, t = line.split("\t")
                    N_ED_MJ += [int(n)]; T_ED_MJ += [float(t)]
                N_ED_MJ, T_ED_MJ = sort_data(N_ED_MJ, T_ED_MJ)
                N_ED_MJ_extrap, T_ED_MJ_extrap = extrapolate_data(N_ED_MJ, T_ED_MJ, "ed mj step " + stepsize)
                # plotting stepsize influence
                subfig1.plot(N_ED_MJ, T_ED_MJ, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_ED_MJ_extrap, T_ED_MJ_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED-fit: Schrittweite %.2f" % float(stepsize))

                # QT data
                QT_file = open("./results/benchmarking/runtime/data/" + "QT_SG_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".txt")
                QT_lines = QT_file.readlines()
                QT_file.close()
                N_QT = []; T_QT = []
                for line in QT_lines:
                    n, t = line.split("\t")
                    N_QT += [int(n)]; T_QT += [float(t)]
                N_QT, T_QT = sort_data(N_QT, T_QT)
                N_QT_extrap, T_QT_extrap = extrapolate_data(N_QT, T_QT, "qt, step " + stepsize)
                # plotting stepsize influence
                subfig1.plot(N_QT, T_QT, lw = 0.0, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_QT_extrap, T_QT_extrap, lw = line_width, ls = "dashed", markersize = 0.0, marker = "o", color = colors[color_count], label = "QT-fit:  Schrittweite %.2f" % float(stepsize))
                color_count += 1
    
            subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = labelfontsize)
            subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = labelfontsize)
            if sample == "1": vec = "Startvektor"
            else: vec = "Startvektoren"
            subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$" + "\n" + "mit %s %s bei der QT" % (sample, vec), fontsize = titlefontsize)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
            subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
            add_time_steps(subfig1)
            fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
            fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".pdf")
            subfig1.set_yscale('log')
            fig1.savefig("./results/benchmarking/runtime/" + "log_QT_ED_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
            fig1.savefig("./results/benchmarking/runtime/" + "log_QT_ED_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".pdf")
            plt.close(fig1)

##### run time, m_z = 0 block #####
def RT_plot_raw_files_mag_zero():
    print("plotting raw files (m_z = 0 block)")
    for core in cores:
        if core != "1": continue
        # ED exac data
        # momentum states
        ED_MS_file = open("./results/benchmarking/runtime/data/" + "ED_MS_" + core + ".txt")
        ED_MS_lines = ED_MS_file.readlines()
        ED_MS_file.close()
        N_ED_MS = []; T_ED_MS = []
        for line in ED_MS_lines:
            n, t = line.split("\t")
            N_ED_MS += [int(n)]; T_ED_MS += [float(t)]
        N_ED_MS, T_ED_MS = sort_data(N_ED_MS, T_ED_MS)
        N_ED_MS_extrap, T_ED_MS_extrap = extrapolate_data(N_ED_MS, T_ED_MS, "ed ms")

        # spin inversion
        ED_SI_file = open("./results/benchmarking/runtime/data/" + "ED_SI_" + core + ".txt")
        ED_SI_lines = ED_SI_file.readlines()
        ED_SI_file.close()
        N_ED_SI = []; T_ED_SI = []
        for line in ED_SI_lines:
            n, t = line.split("\t")
            N_ED_SI += [int(n)]; T_ED_SI += [float(t)]
        N_ED_SI, T_ED_SI = sort_data(N_ED_SI, T_ED_SI)
        N_ED_SI_extrap, T_ED_SI_extrap = extrapolate_data(N_ED_SI, T_ED_SI, "ed si")

        for stepsize in stepsizes:

            for sample in samples:
                # QT data
                QT_file = open("./results/benchmarking/runtime/data/" + "QT_SG_zero_block_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".txt")
                QT_lines = QT_file.readlines()
                QT_file.close()
                N_QT = []; T_QT = []
                for line in QT_lines:
                    n, t = line.split("\t")
                    N_QT += [int(n)]; T_QT += [float(t)]
                N_QT, T_QT = sort_data(N_QT, T_QT)
                N_QT_extrap, T_QT_extrap = extrapolate_data(N_QT, T_QT, "qt")

                # plotting
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                color_count = 0
                subfig1.plot(N_ED_MS, T_ED_MS, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_ED_MS_extrap, T_ED_MS_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED: MS")
                color_count += 1
                subfig1.plot(N_ED_SI, T_ED_SI, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_ED_SI_extrap, T_ED_SI_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED: SI")
                color_count += 1
                subfig1.plot(N_QT, T_QT, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
                subfig1.plot(N_QT_extrap, T_QT_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "QT: MS")
                color_count += 1
                subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
                subfig1.set_ylabel(r'$t$ in $s$', fontsize = labelfontsize)
                if sample == "1": vec = "Startvektor"
                else: vec = "Startvektoren"
                subfig1.set_title(r"Laufzeit $t$ des $m_z = 0$ Blocks für unterschiedliche Systemgrößen $N$" + "\n" + "nur QT: Schrittweite = %.2f, %s %s" % (float(stepsize), sample, vec), fontsize = titlefontsize)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
                subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
                add_time_steps(subfig1)
                fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_SG_MS_SI_QT_MS_zero_block_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
                fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_SG_MS_SI_QT_MS_zero_block_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".pdf")
                subfig1.set_yscale('log')
                fig1.savefig("./results/benchmarking/runtime/" + "log_QT_ED_SG_MS_SI_QT_MS_zero_block_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
                fig1.savefig("./results/benchmarking/runtime/" + "log_QT_ED_SG_MS_SI_QT_MS_zero_block_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".pdf")
                plt.close(fig1)

def RT_plot_only_ED_mag_zero():
    print("comparing only ED (m_z = 0 block)")
    for core in cores:
        # momentum states
        ED_MS_file = open("./results/benchmarking/runtime/data/" + "ED_MS_" + core + ".txt")
        ED_MS_lines = ED_MS_file.readlines()
        ED_MS_file.close()
        N_ED_MS = []; T_ED_MS = []
        for line in ED_MS_lines:
            n, t = line.split("\t")
            N_ED_MS += [int(n)]; T_ED_MS += [float(t)]
        N_ED_MS, T_ED_MS = sort_data(N_ED_MS, T_ED_MS)
        N_ED_MS_extrap, T_ED_MS_extrap = extrapolate_data(N_ED_MS, T_ED_MS, "ed only ms")

        # spin inversion
        ED_SI_file = open("./results/benchmarking/runtime/data/" + "ED_SI_" + core + ".txt")
        ED_SI_lines = ED_SI_file.readlines()
        ED_SI_file.close()
        N_ED_SI = []; T_ED_SI = []
        for line in ED_SI_lines:
            n, t = line.split("\t")
            N_ED_SI += [int(n)]; T_ED_SI += [float(t)]
        N_ED_SI, T_ED_SI = sort_data(N_ED_SI, T_ED_SI)
        N_ED_SI_extrap, T_ED_SI_extrap = extrapolate_data(N_ED_SI, T_ED_SI,  "ed only si")

        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        color_count = 0
        subfig1.plot(N_ED_MS, T_ED_MS, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
        subfig1.plot(N_ED_MS_extrap, T_ED_MS_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED: MS")
        color_count += 1
        subfig1.plot(N_ED_SI, T_ED_SI, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
        subfig1.plot(N_ED_SI_extrap, T_ED_SI_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = "ED: SI")
        color_count += 1

        subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = labelfontsize)
        subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = labelfontsize)
        subfig1.set_title(r"Laufzeit $t$ des $m_z = 0$ Blocks für unterschiedliche Systemgrößen $N$", fontsize = titlefontsize)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        add_time_steps(subfig1)
        fig1.savefig("./results/benchmarking/runtime/" + "ED_SG_MS_SI_zero_block" + "_cores_" + core + ".png")
        fig1.savefig("./results/benchmarking/runtime/" + "ED_SG_MS_SI_zero_block" + "_cores_" + core + ".pdf")
        subfig1.set_yscale('log')
        fig1.savefig("./results/benchmarking/runtime/" + "log_ED_SG_MS_SI_zero_block" + "_cores_" + core + ".png")
        fig1.savefig("./results/benchmarking/runtime/" + "log_ED_SG_MS_SI_zero_block" + "_cores_" + core + ".pdf")
        plt.close(fig1)

##### memory usage #####
def MU_plot_raw_files():
    print("plotting raw files")
    # QT H
    file = open("results/benchmarking/memoryusage/data/QT_H.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "QT H")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = r"QT")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$ bei der QT", fontsize = titlefontsize)
    subfig1.set_xticks(N)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig1.savefig("./results/benchmarking/memoryusage/QT_H.png")
    fig1.savefig("./results/benchmarking/memoryusage/QT_H.pdf")
    subfig1.set_yscale("log")
    fig1.savefig("./results/benchmarking/memoryusage/log_QT_H.png")
    fig1.savefig("./results/benchmarking/memoryusage/log_QT_H.pdf")
    plt.close(fig1)

    # QT S2
    file = open("results/benchmarking/memoryusage/data/QT_S2.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "QT S2")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = r"QT")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Spinmatrix $S^2$ bei der QT", fontsize = titlefontsize)
    subfig1.set_xticks(N)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig1.savefig("./results/benchmarking/memoryusage/QT_S2.png")
    fig1.savefig("./results/benchmarking/memoryusage/QT_S2.pdf")
    subfig1.set_yscale("log")
    fig1.savefig("./results/benchmarking/memoryusage/log_QT_S2.png")
    fig1.savefig("./results/benchmarking/memoryusage/log_QT_S2.pdf")
    plt.close(fig1)

    # ED H S2 MS
    file = open("results/benchmarking/memoryusage/data/ED_H_S2_MS.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED H S2")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = r"ED")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamilton- $H$ \& Spinmatrix $S^2$ bei der ED (Impulszutände)", fontsize = titlefontsize)
    subfig1.set_xticks(N)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig1.savefig("./results/benchmarking/memoryusage/ED_H_S2_MS.png")
    fig1.savefig("./results/benchmarking/memoryusage/ED_H_S2_MS.pdf")
    subfig1.set_yscale("log")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_H_S2_MS.png")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_H_S2_MS.pdf")
    plt.close(fig1)

    # ED H SI
    file = open("results/benchmarking/memoryusage/data/ED_H_SI.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED H SI")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = r"ED")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$ bei der ED (Spininversion)", fontsize = titlefontsize)
    subfig1.set_xticks(N)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig1.savefig("./results/benchmarking/memoryusage/ED_H_SI.png")
    fig1.savefig("./results/benchmarking/memoryusage/ED_H_SI.pdf")
    subfig1.set_yscale("log")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_H_SI.png")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_H_SI.pdf")
    plt.close(fig1)

def MU_plot_only_ED():
    print("comparing only ED")

    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

    # ED H S2 MS
    file = open("results/benchmarking/memoryusage/data/ED_H_S2_MS.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED H S2 comp")
    color_count = 0
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = r"MS")
    color_count += 1

    # ED H SI
    file = open("results/benchmarking/memoryusage/data/ED_H_SI.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED H SI comp")
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count], label = r"SI")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$, Spininversion (SI) \& Impulszustände (MS)", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig1.savefig("./results/benchmarking/memoryusage/ED_MS_SI.png")
    fig1.savefig("./results/benchmarking/memoryusage/ED_MS_SI.pdf")
    subfig1.set_yscale("log")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_MS_SI.png")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_MS_SI.pdf")
    plt.close(fig1)

def MU_ED_vs_QT():
    print("comparing ED and QT")

    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count1 = 0
    fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
    color_count2 = 0
    fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
    color_count3 = 0

    x_start = 5; x_end = 33
    x_ticks = [6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]

    ##### get data & plot #####

    # ED H S2 MS
    file = open("results/benchmarking/memoryusage/data/ED_H_S2_MS.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED QT comp qt S2 ED MS")
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count1])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count1], label = r"\textit{dense} (MS)")
    color_count1 += 1
    subfig2.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count2])
    subfig2.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count2], label = r"\textit{dense} (MS)")
    color_count2 += 1
    subfig3.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count3])
    subfig3.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count3], label = r"\textit{dense} (MS): $H$ \& $S^2$")
    color_count3 += 1

    # ED H S2 SI
    file = open("results/benchmarking/memoryusage/data/ED_H_SI.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)/2.0]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED QT comp qt S2 ED SI")
    subfig1.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count1])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count1], label = r"\textit{dense} (SI)$^*$")
    color_count1 += 1
    subfig2.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count2])
    subfig2.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count2], label = r"\textit{dense} (SI)$^*$")
    color_count2 += 1
    subfig3.plot(N, RAM, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count3])
    subfig3.plot(N_extrap, RAM_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count3], label = r"\textit{dense} (SI)$^*$: $H$ \& $S^2$")
    color_count3 += 1

    # QT H
    file = open("results/benchmarking/memoryusage/data/QT_H.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED QT comp qt H")
    subfig1.plot(N, RAM, lw = 0.0, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count1])
    subfig1.plot(N_extrap, RAM_extrap, lw = line_width, ls = "dashed", markersize = 0.0, marker = "o", color = colors[color_count1], label = r"\textit{sparse} (MS)")
    color_count1 += 1
    subfig3.plot(N, RAM, lw = 0.0, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count3])
    subfig3.plot(N_extrap, RAM_extrap, lw = line_width, ls = "dashed", markersize = 0.0, marker = "o", color = colors[color_count3], label = r"\textit{sparse} (MS): $H$")
    color_count3 += 1

    # QT S2
    file = open("results/benchmarking/memoryusage/data/QT_S2.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    N_extrap, RAM_extrap = extrapolate_data(N, RAM, "ED QT comp qt S2")
    subfig2.plot(N, RAM, lw = 0.0, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count2])
    subfig2.plot(N_extrap, RAM_extrap, lw = line_width, ls = "dashed", markersize = 0.0, marker = "o", color = colors[color_count2], label = r"\textit{sparse} (MS)")
    color_count2 += 1
    subfig3.plot(N, RAM, lw = 0.0, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count3])
    subfig3.plot(N_extrap, RAM_extrap, lw = line_width, ls = "dashed", markersize = 0.0, marker = "o", color = colors[color_count3], label = r"\textit{sparse} (MS): $S^2$")
    color_count3 += 1

    ##### save plots #####

    add_size_steps(subfig1)
    add_size_steps(subfig2)
    add_size_steps(subfig3)

    subfig1.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$, \textit{dense}- \& \textit{sparse}-Matrizen", fontsize = titlefontsize)
    subfig1.set_xticks(N)
    subfig1.set_xlim(x_start, x_end)
    subfig1.set_xticks(x_ticks)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig1.savefig("./results/benchmarking/memoryusage/ED_vs_QT_H.png")
    fig1.savefig("./results/benchmarking/memoryusage/ED_vs_QT_H.pdf")
    subfig1.set_yscale('log', base = 2)
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_vs_QT_H.png")
    fig1.savefig("./results/benchmarking/memoryusage/log_ED_vs_QT_H.pdf")
    plt.close(fig1)

    subfig2.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig2.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig2.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Spinmatrix $S^2$, \textit{dense}- \& \textit{sparse}-Matrizen", fontsize = titlefontsize)
    subfig2.set_xticks(N)
    subfig2.set_xlim(x_start, x_end)
    subfig2.set_xticks(x_ticks)
    subfig2.axhline(0, color = "grey")
    subfig2.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
    subfig2.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig2.savefig("./results/benchmarking/memoryusage/ED_vs_QT_S2.png")
    fig2.savefig("./results/benchmarking/memoryusage/ED_vs_QT_S2.pdf")
    subfig2.set_yscale('log', base = 2)
    fig2.savefig("./results/benchmarking/memoryusage/log_ED_vs_QT_S2.png")
    fig2.savefig("./results/benchmarking/memoryusage/log_ED_vs_QT_S2.pdf")
    plt.close(fig2)

    subfig3.set_xlabel(r'$N$', fontsize = labelfontsize)
    subfig3.set_ylabel(r'\# Matrixelemente', fontsize = labelfontsize)
    subfig3.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamilton- $H$ \& Spinmatrix $S^2$, \textit{dense}- \& \textit{sparse}-Matrizen", fontsize = titlefontsize)
    subfig3.set_xticks(N)
    subfig3.set_xlim(x_start, x_end)
    subfig3.set_xticks(x_ticks)
    subfig3.axhline(0, color = "grey")
    subfig3.legend(loc = 'lower right' ,frameon = False, fontsize = legendfontsize)
    subfig3.tick_params(axis="both", which="major", labelsize=axisfontsize)
    fig3.savefig("./results/benchmarking/memoryusage/ED_vs_QT_H_S2.png")
    fig3.savefig("./results/benchmarking/memoryusage/ED_vs_QT_H_S2.pdf")
    subfig3.set_yscale('log', base =2 )
    fig3.savefig("./results/benchmarking/memoryusage/log_ED_vs_QT_H_S2.png")
    fig3.savefig("./results/benchmarking/memoryusage/log_ED_vs_QT_H_S2.pdf")
    plt.close(fig3)

##### main #####
if __name__ == "__main__":

    # print("plotting bechmarking (run time):")
    # RT_plot_raw_files()
    # RT_plot_only_ED()
    # RT_plot_step_size_influence()
    RT_plot_raw_files_mag_zero()
    RT_plot_only_ED_mag_zero()

    print("plotting bechmarking (memory usage):")
    MU_plot_raw_files()
    MU_plot_only_ED()
    MU_ED_vs_QT()

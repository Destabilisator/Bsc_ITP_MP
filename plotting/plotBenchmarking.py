#-*- coding:utf-8 -*-
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

colors = ["red", "blue", "green", "tomato", "purple"]

line_width = 3
marker_size = 5

counter = 0

samples = ["1"]#, "2", "3"]
stepsizes = ["0.100000", "0.010000"]
cores = ["1"] # ["1", "2", "5", "10"]

def sort_data(N, T):
    length = len(N)
    for i in range(length-1):
        for j in range(0, length-i-1):
            N1 = N[j]; N2 = N[j+1]
            if N1 > N2:
                N[j], N[j+1] = N[j+1], N[j]
                T[j], T[j+1] = T[j+1], T[j]
    return N, T

def extrap_func(x: float, A: float, k: float) -> float:
#def extrap_func(x, A, k):
    return A * np.exp(k * x)

def extrapolate_data(N, T):
    fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
    subfig2.plot(N, T, lw = 1, ls = "solid", markersize = 5, marker = "o", color = "black")
    fig2.savefig("./results/benchmarking/temp/" + str(counter) + ".png")
    counter += 1

    params, cv = scipy.optimize.curve_fit(extrap_func, N, T, (0.1, 0.1))
    A_param, k_param = params
    X = np.linspace(6, 32, 1000)
    Y = extrap_func(X, A_param, k_param)
    return X, Y

# run time, fll matrix
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

                # plotting
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                color_count = 0
                subfig1.plot(N_ED_SG, T_ED_SG, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-ev")
                color_count += 1
                subfig1.plot(N_ED_MJ, T_ED_MJ, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-fit")
                color_count += 1
                subfig1.plot(N_QT, T_QT, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "QT-fit")
                color_count += 1
                subfig1.set_xlabel(r'$N$', fontsize = 40)
                subfig1.set_ylabel(r'$t$ in $s$', fontsize = 40)
                if sample == "1": vec = "Startvektor"
                else: vec = "Startvektoren"
                subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$" + "\n" + "Schrittweite (nur Fits) = %.2f, %s %s (nur QT)" % (float(stepsize), sample, vec), fontsize = 40)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
                fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_SG_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
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

        subfig1.plot(N_ED_SG, T_ED_SG, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-ev")
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
            N_ED_MJ, T_ED_MJ = sort_data(N_ED_MJ, T_ED_MJ)

            subfig1.plot(N_ED_MJ, T_ED_MJ, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-fit: Schrittweite %.2f" % float(stepsize))
            color_count += 1

        subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = 40)
        subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = 40)
        subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$", fontsize = 40)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
        fig1.savefig("./results/benchmarking/runtime/" + "ED_SG_MJ" + "_cores_" + core + ".png")
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
                # plotting stepsize influence
                subfig1.plot(N_ED_MJ, T_ED_MJ, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-fit: Schrittweite %.2f" % float(stepsize))

                # QT data
                QT_file = open("./results/benchmarking/runtime/data/" + "QT_SG_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".txt")
                QT_lines = QT_file.readlines()
                QT_file.close()
                N_QT = []; T_QT = []
                for line in QT_lines:
                    n, t = line.split("\t")
                    N_QT += [int(n)]; T_QT += [float(t)]
                N_QT, T_QT = sort_data(N_QT, T_QT)
                # plotting stepsize influence
                subfig1.plot(N_QT, T_QT, lw = line_width, ls = "dashed", markersize = marker_size, marker = "o", color = colors[color_count], label = "QT-fit:  Schrittweite %.2f" % float(stepsize))
                color_count += 1
    
            subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = 40)
            subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = 40)
            if sample == "1": vec = "Startvektor"
            else: vec = "Startvektoren"
            subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$" + "\n" + "mit %s %s bei der QT" % (sample, vec), fontsize = 40)
            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
            fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
            plt.close(fig1)

# run time, m_z = 0 block
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
        N_ED_MS_extrap, T_ED_MS_extrap = extrapolate_data(N_ED_MS, T_ED_MS)

        # spin inversion
        ED_SI_file = open("./results/benchmarking/runtime/data/" + "ED_SI_" + core + ".txt")
        ED_SI_lines = ED_SI_file.readlines()
        ED_SI_file.close()
        N_ED_SI = []; T_ED_SI = []
        for line in ED_SI_lines:
            n, t = line.split("\t")
            N_ED_SI += [int(n)]; T_ED_SI += [float(t)]
        N_ED_SI, T_ED_SI = sort_data(N_ED_SI, T_ED_SI)
        N_ED_SI_extrap, T_ED_SI_extrap = extrapolate_data(N_ED_SI, T_ED_SI)

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
                N_QT_extrap, T_QT_extrap = extrapolate_data(N_QT, T_QT)

                # plotting
                fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
                color_count = 0
                subfig1.plot(N_ED_MS, T_ED_MS, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED: MS")
                subfig1.plot(N_ED_MS_extrap, T_ED_MS_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count])
                color_count += 1
                subfig1.plot(N_ED_SI, T_ED_SI, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED: SI")
                subfig1.plot(N_ED_SI_extrap, T_ED_SI_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count])
                color_count += 1
                subfig1.plot(N_QT, T_QT, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "QT: MS")
                subfig1.plot(N_QT, T_QT_extrap, lw = line_width_extrap, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count])
                color_count += 1
                subfig1.set_xlabel(r'$N$', fontsize = 40)
                subfig1.set_ylabel(r'$t$ in $s$', fontsize = 40)
                if sample == "1": vec = "Startvektor"
                else: vec = "Startvektoren"
                subfig1.set_title(r"Laufzeit $t$ des $m_z = 0$ Blocks für unterschiedliche Systemgrößen $N$" + "\n" + "nur QT: Schrittweite = %.2f, %s %s" % (float(stepsize), sample, vec), fontsize = 40)
                subfig1.axhline(0, color = "grey")
                subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
                fig1.savefig("./results/benchmarking/runtime/" + "QT_ED_SG_MS_SI_QT_MS_zero_block_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
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
        N_ED_MS_extrap, T_ED_MS_extrap = extrapolate_data(N_ED_MS, T_ED_MS)

        # spin inversion
        ED_SI_file = open("./results/benchmarking/runtime/data/" + "ED_SI_" + core + ".txt")
        ED_SI_lines = ED_SI_file.readlines()
        ED_SI_file.close()
        N_ED_SI = []; T_ED_SI = []
        for line in ED_SI_lines:
            n, t = line.split("\t")
            N_ED_SI += [int(n)]; T_ED_SI += [float(t)]
        N_ED_SI, T_ED_SI = sort_data(N_ED_SI, T_ED_SI)
        N_ED_SI_extrap, T_ED_SI_extrap = extrapolate_data(N_ED_SI, T_ED_SI)

        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        color_count = 0
        subfig1.plot(N_ED_MS, T_ED_MS, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED: MS")
        subfig1.plot(N_ED_MS_extrap, T_ED_MS_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count])
        color_count += 1
        subfig1.plot(N_ED_SI, T_ED_SI, lw = 0.0, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED: SI")
        subfig1.plot(N_ED_SI_extrap, T_ED_SI_extrap, lw = line_width, ls = "solid", markersize = 0.0, marker = "o", color = colors[color_count])
        color_count += 1

        subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = 40)
        subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = 40)
        subfig1.set_title(r"Laufzeit $t$ des $m_z = 0$ Blocks für unterschiedliche Systemgrößen $N$", fontsize = 40)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
        fig1.savefig("./results/benchmarking/runtime/" + "ED_SG_MS_SI_zero_block" + "_cores_" + core + ".png")
        plt.close(fig1)

# memory usage
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
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = r"QT")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$ bei der QT", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.set_yscale("log")
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig1.savefig("./results/benchmarking/memoryusage/QT_H.png")
    plt.close(fig1)

    # QT S2
    file = open("results/benchmarking/memoryusage/data/QT_S2.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = r"QT")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Spinmatrix $S^2$ bei der QT", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.set_yscale("log")
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig1.savefig("./results/benchmarking/memoryusage/QT_S2.png")
    plt.close(fig1)

    # ED H S2 MS
    file = open("results/benchmarking/memoryusage/data/ED_H_S2_MS.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = r"ED")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamilton- $H$ \& Spinmatrix $S^2$ bei der ED (Impulszutände)", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.set_yscale("log")
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig1.savefig("./results/benchmarking/memoryusage/ED_H_S2_MS.png")
    plt.close(fig1)

    # ED H SI
    file = open("results/benchmarking/memoryusage/data/ED_H_SI.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count = 0
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = r"ED")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$ bei der ED (Spininversion)", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.set_yscale("log")
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig1.savefig("./results/benchmarking/memoryusage/ED_H_SI.png")
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
    color_count = 0
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = r"MS")
    color_count += 1

    # ED H SI
    file = open("results/benchmarking/memoryusage/data/ED_H_SI.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = r"SI")
    color_count += 1
    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$, Spininversion (SI) \& Impulszustände (MS)", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.set_yscale("log")
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig1.savefig("./results/benchmarking/memoryusage/ED_MS_SI.png")
    plt.close(fig1)

def MU_ED_vs_QT():
    print("comparing ED and QT")

    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    color_count1 = 0
    fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
    color_count2 = 0
    fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
    color_count3 = 0

    # QT H
    file = open("results/benchmarking/memoryusage/data/QT_H.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count1], label = r"QT")
    color_count1 += 1
    subfig3.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count3], label = r"QT: $H$")
    color_count3 += 1

    # QT S2
    file = open("results/benchmarking/memoryusage/data/QT_S2.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    subfig2.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count2], label = r"QT")
    color_count2 += 1
    subfig3.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count3], label = r"QT: $S^2$")
    color_count3 += 1

    # ED H S2 MS
    file = open("results/benchmarking/memoryusage/data/ED_H_S2_MS.txt")
    lines = file.readlines()
    file.close()
    N = []; RAM = []
    for line in lines:
        n, ram = line.split("\t")
        N += [int(n)]; RAM += [float(ram)]
    subfig1.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count1], label = r"ED")
    color_count1 += 1
    subfig2.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count2], label = r"ED")
    color_count2 += 1
    subfig3.plot(N, RAM, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count3], label = r"ED: $H$ \& $S^2$")
    color_count3 += 1

    subfig1.set_xlabel(r'$N$', fontsize = 40)
    subfig1.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig1.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamiltonmatrix $H$, QT \& ED (Impulszutände)", fontsize = 40)
    subfig1.set_xticks(N)
    subfig1.set_yscale("log")
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig1.savefig("./results/benchmarking/memoryusage/ED_vs_QT_H.png")
    plt.close(fig1)

    subfig2.set_xlabel(r'$N$', fontsize = 40)
    subfig2.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig2.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Spinmatrix $S^2$, QT \& ED (Impulszutände)", fontsize = 40)
    subfig2.set_xticks(N)
    subfig2.set_yscale("log")
    subfig2.axhline(0, color = "grey")
    subfig2.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig2.savefig("./results/benchmarking/memoryusage/ED_vs_QT_S2.png")
    plt.close(fig2)

    subfig3.set_xlabel(r'$N$', fontsize = 40)
    subfig3.set_ylabel(r'\# Matrixelemente', fontsize = 40)
    subfig3.set_title(r"Anzahl der zu speichernden Matrixelemente" + "\n" + r"Hamilton- $H$ \& Spinmatrix $S^2$, QT \& ED (Impulszutände)", fontsize = 40)
    subfig3.set_xticks(N)
    subfig3.set_yscale("log")
    subfig3.axhline(0, color = "grey")
    subfig3.legend(loc = 'best' ,frameon = False, fontsize = 30)
    fig3.savefig("./results/benchmarking/memoryusage/ED_vs_QT_H_S2.png")
    plt.close(fig3)


if __name__ == "__main__":

    print("plotting bechmarking (run time):")
    RT_plot_raw_files()
    RT_plot_only_ED()
    RT_plot_step_size_influence()
    RT_plot_raw_files_mag_zero()
    RT_plot_only_ED_mag_zero()

    print("plotting bechmarking (memory usage):")
    MU_plot_raw_files()
    MU_plot_only_ED()
    MU_ED_vs_QT()

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
plt.rcParams['text.usetex'] = True

colors = ["red", "blue", "green", "tomato", "purple"]

line_width = 3
marker_size = 5

samples = ["1", "2", "3"]
stepsizes = ["0.100000", "0.010000"]
cores = ["1"] # ["1", "2", "5", "10"]

def plot_raw_files():
    print("plotting raw files")
    for core in cores:
        # ED exac data
        ED_SG_file = open("./results/benchmarking/data/" + "ED_SG_cores_" + core + ".txt")
        ED_SG_lines = ED_SG_file.readlines()
        ED_SG_file.close()
        N_ED_SG = []; T_ED_SG = []
        for line in ED_SG_lines:
            n, t = line.split("\t")
            N_ED_SG += [int(n)]; T_ED_SG += [float(t)]

        for stepsize in stepsizes:
            # ED fit data
            ED_MJ_file = open("./results/benchmarking/data/" + "ED_MJ_step_" + stepsize + "_cores_" + core + ".txt")
            ED_MJ_lines = ED_MJ_file.readlines()
            ED_MJ_file.close()
            N_ED_MJ = []; T_ED_MJ = []
            for line in ED_MJ_lines:
                n, t = line.split("\t")
                N_ED_MJ += [int(n)]; T_ED_MJ += [float(t)]

            for sample in samples:
                # QT data
                QT_file = open("./results/benchmarking/data/" + "QT_SG_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".txt")
                QT_lines = QT_file.readlines()
                QT_file.close()
                N_QT = []; T_QT = []
                for line in QT_lines:
                    n, t = line.split("\t")
                    N_QT += [int(n)]; T_QT += [float(t)]

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
                fig1.savefig("./results/benchmarking/" + "QT_ED_SG_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
                plt.close(fig1)

def plot_only_ED():
    print("comparing only ED")
    for core in cores:
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        color_count = 0
        # ED exac data
        ED_SG_file = open("./results/benchmarking/data/" + "ED_SG_cores_" + core + ".txt")
        ED_SG_lines = ED_SG_file.readlines()
        ED_SG_file.close()
        N_ED_SG = []; T_ED_SG = []
        for line in ED_SG_lines:
            n, t = line.split("\t")
            N_ED_SG += [int(n)]; T_ED_SG += [float(t)]

        subfig1.plot(N_ED_SG, T_ED_SG, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-ev")
        color_count += 1

        for stepsize in stepsizes:
            # ED fit data
            ED_MJ_file = open("./results/benchmarking/data/" + "ED_MJ_step_" + stepsize + "_cores_" + core + ".txt")
            ED_MJ_lines = ED_MJ_file.readlines()
            ED_MJ_file.close()
            N_ED_MJ = []; T_ED_MJ = []
            for line in ED_MJ_lines:
                n, t = line.split("\t")
                N_ED_MJ += [int(n)]; T_ED_MJ += [float(t)]

            subfig1.plot(N_ED_MJ, T_ED_MJ, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-fit: Schrittweite %.2f" % float(stepsize))
            color_count += 1

        subfig1.set_xlabel(r'Gitterplätze $N$', fontsize = 40)
        subfig1.set_ylabel(r'$t$ in Sekunden', fontsize = 40)
        subfig1.set_title(r"Laufzeit $t$ für unterschiedliche Systemgrößen $N$", fontsize = 40)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)
        fig1.savefig("./results/benchmarking/" + "ED_SG_MJ" + "_cores_" + core + ".png")
        plt.close(fig1)

def plot_step_size_influence():
    print("comparing influence of stepsize")
    for core in cores:
        for sample in samples:
            # fig for stepsize influence
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            color_count = 0
            for stepsize in stepsizes:
                # ED fit data
                ED_MJ_file = open("./results/benchmarking/data/" + "ED_MJ_step_" + stepsize + "_cores_" + core + ".txt")
                ED_MJ_lines = ED_MJ_file.readlines()
                ED_MJ_file.close()
                N_ED_MJ = []; T_ED_MJ = []
                for line in ED_MJ_lines:
                    n, t = line.split("\t")
                    N_ED_MJ += [int(n)]; T_ED_MJ += [float(t)]
                # plotting stepsize influence
                subfig1.plot(N_ED_MJ, T_ED_MJ, lw = line_width, ls = "solid", markersize = marker_size, marker = "o", color = colors[color_count], label = "ED-fit: Schrittweite %.2f" % float(stepsize))

                # QT data
                QT_file = open("./results/benchmarking/data/" + "QT_SG_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".txt")
                QT_lines = QT_file.readlines()
                QT_file.close()
                N_QT = []; T_QT = []
                for line in QT_lines:
                    n, t = line.split("\t")
                    N_QT += [int(n)]; T_QT += [float(t)]
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
            fig1.savefig("./results/benchmarking/" + "QT_ED_MJ_step_" + stepsize + "_SAMPLES_" + sample + "_cores_" + core + ".png")
            plt.close(fig1)


if __name__ == "__main__":
    print("plotting bechmarking:")

    plot_raw_files()
    plot_only_ED()
    plot_step_size_influence()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import gc
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color = [("12", "magenta"), ("14", "brown"), ("16", "purple")]#, ("18", "tomato")]
N_color = [("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]

colors = ["red", "blue", "green", "magenta", "brown", "purple", "tomato"]

max_n = 30

einheit_x = r'$k_B T$ / $J_2$' #'$T$ in $k_B$ / $J_2$'

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

start = 0.0
ends = [0.1, 0.15, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0]#, 20.0, 50.0, 100.0]

print("plotting susceptibility (constant J1/J2, funtion of T) ...")
N_color = [("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color = [("16", "purple")]
for N_outer in range(len(N_color)):
    continue
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
    # X_high_T = np.linspace(0.01, end, 5000)
    # Y_high_T = 3/4 / X_high_T**2
    used_N = "N"
    for N_inner in range(N_outer, len(N_color)):
        x_min_ED = 42069.0; x_min_QT = 42069.0
        print("N_outer: %s, N_inner: %s" % (N_color[N_outer][0], N_color[N_inner][0]))
        fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
        # ED
        N, c = N_color [N_inner]
        file = open("results/" + N + "/data/data_susceptibility_J_const.txt", 'r')
        lines = file.readlines()
        linesJ = lines[0][len("J1/J2: "):-1]
        lbl = "N = " + N
        X = []
        Y = []
        for i in range(8,len(lines)):
            x, y = lines[i].split("\t")
            # if float(x) <= start or float(x) > end: continue
            if float(x) <= 0.0: continue
            X += [1/float(x)]
            Y += [float(y)]
        file.close()
        if min(X) < x_min_ED: x_min_ED = min(X)
        subfig1.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = lbl)
        subfig3.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = "ED")

        # QT
        file = open("results/" + N + "/data/1_data_susceptibility_J_const_QT.txt", 'r')
        lines = file.readlines()
        linesJ = lines[1][len("J1/J2: "):-1]
        lbl = "N = " + N
        X = []
        Y = []
        YErr = []
        for i in range(6,len(lines)):
            x, y, yErr = lines[i].split("\t")
            # if float(x) <= start or float(x) > end : continue
            if float(x) <= 0.0: continue
            X += [1/float(x)]
            Y += [float(y)]
            YErr += [float(yErr)]
        file.close()
        if min(X) < x_min_QT: x_min_QT = min(X)
        subfig3.plot(X, Y, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = c, label = "QT")
        X = np.asarray(X)
        Y = np.asarray(Y)
        YErr = np.asarray(YErr)
        subfig3.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.20)

        # subfig3.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 40)
        subfig3.set_xlabel(einheit_x, fontsize = labelfontsize)
        subfig3.set_ylabel('$\\chi/N$ / $J_2$', fontsize = labelfontsize)
        subfig3.set_title("Suszeptibilität pro Spin $\\chi/N$ bei N = " + N + "\nbei 12 Mittelungen über einen Startvektor", fontsize = titlefontsize)

        subfig3.axhline(0, color = "grey")
        subfig3.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

        # plt.xticks(fontsize = axisfontsize)
        # plt.yticks(fontsize = axisfontsize)
        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)

        # subfig3.set_ylim(0.0, 0.5)
        # subfig3.set_xscale("log")

        for end in ends:
            subfig3.set_xlim(min(x_min_ED, x_min_QT), end)
            fig3.savefig("results/" + "X_ED_QT_" + N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            fig3.savefig("results/" + "X_ED_QT_" + N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        plt.close(fig3)

        subfig2.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = lbl)
        subfig2.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.15)

        used_N += "_" + N

    # subfig1.set_xlabel(r'$\beta$ in $k_B$ / $J_2$', fontsize = 40)
    subfig1.set_xlabel(einheit_x, fontsize = labelfontsize)
    subfig1.set_ylabel('$\\chi/N$ / $J_2$', fontsize = labelfontsize)
    subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$', fontsize = titlefontsize)

    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

    # plt.xticks(fontsize = 25)
    # plt.yticks(fontsize = 25)
    subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
    #subfig1.tick_params(axis="both", which="major", labelsize=25)

    # subfig1.set_ylim(0.0, 0.5)
    # subfig1.set_xscale("log")

    for end in ends:
        subfig1.set_xlim(x_min_ED, end)
        fig1.savefig("results/" + "X_ED_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
        fig1.savefig("results/" + "X_ED_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")

    # subfig1.plot(X_high_T, Y_high_T, lw = 4, ls = "solid", markersize = 0, marker = "o", color = "black", label = "high T")
    # fig1.savefig("results/" + "X_ED_high_T_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")

    plt.close(fig1)


    #subfig2.set_xlabel(r'$\beta$ in $k_B$ / $J_2$', fontsize = 40)
    subfig2.set_xlabel(einheit_x, fontsize = labelfontsize)
    subfig2.set_ylabel('$\\chi/N$ / $J_2$', fontsize = labelfontsize)
    subfig2.set_title('Suszeptibilität pro Spin $\\chi/N$' + "\nmit 12 Mittelungen über einen Startvektor", fontsize = titlefontsize)

    subfig2.axhline(0, color = "grey")
    subfig2.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

    # plt.xticks(fontsize = 25)
    # plt.yticks(fontsize = 25)
    subfig2.tick_params(axis="both", which="major", labelsize=axisfontsize)

    # subfig2.set_ylim(0.0, 0.5)
    # subfig2.set_xscale("log")

    for end in ends:
        subfig2.set_xlim(x_min_QT, end)
        fig2.savefig("results/" + "X_QT_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
        fig2.savefig("results/" + "X_QT_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")

    # subfig2.plot(X_high_T, Y_high_T, lw = 4, ls = "solid", markersize = 0, marker = "o", color = "black", label = "high T")
    # fig2.savefig("results/" + "X_QT_hight_T_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")

    plt.close(fig2)

print("multiple runs (QT)")
N_color = [("16", "purple"), ("18", "tomato")]
for N_outer in range(len(N_color)):
    N, C = N_color[N_outer]
    for filename in os.listdir("results/" + N + "/data/spin_gap_data/1/"):
        if ".png" in filename or ".pdf" in filename: continue
        if "ED" in filename: continue
        if "placeholder" in filename: continue
        figMultiQT, subfigMultiQT = plt.subplots(1,1,figsize=(16,9))
        J = filename[len("X_J"):-len("QT.txt")]
        print("N: %s, J: %s" % (N_color[N_outer][0], str(J)))
        x_min = 42069
        legend = []
        try:
            filenameED = filename[:-len("QT.txt")] + "ED.txt"
            file = open("results/" + N + "/data/spin_gap_data/1/" + filenameED, 'r')
            lines = file.readlines()
            X = []
            Y = []
            YErr = []
            for i in range(5,len(lines)):
                x, y = lines[i].split("\t")
                #if float(x) < start or float(x) > end: continue
                if float(x) <= 0.0: continue
                X += [1/float(x)]
                Y += [float(y)]
            file.close()
            subfigMultiQT.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = "black")
            if min(X) < x_min: x_min = min(X)
            legend += [Line2D([0], [0], label = "ED", color = "black", ls = "solid", lw = 4)]
        except:
            print("could not plot multiple runs (ED ref) N = %s, J = %s" %(N, J))
        try:
            X_arr = [[]] * max_n; Y_arr = [[]] * max_n
            for n in range(1, max_n+1):
                file = open("results/" + N + "/data/spin_gap_data/" + str(n) + "/" + filename, 'r')
                lines = file.readlines()
                X = []
                Y = []
                YErr = []
                X_arr[n-1] = []; Y_arr[n-1] = []
                for i in range(5,len(lines)):
                    x, y= lines[i].split("\t")
                    #if float(x) < start or float(x) > end: continue
                    if float(x) <= 0.0: continue
                    X += [1/float(x)]; X_arr[n-1] += [1/float(x)]
                    Y += [float(y)]; Y_arr[n-1] += [float(y)]
                file.close()
                subfigMultiQT.plot(X, Y, lw = 2, ls = "solid", markersize = 0, marker = "o", color = "blue", alpha = 0.75)
                if min(X) < x_min: x_min = min(X)
            
            X = []; Y = []; YErr = []
            for pos in range(len(X_arr[0])):
                y = []
                for n in range(max_n):
                    if len(X_arr[n]) == 0: continue
                    y += [Y_arr[n][pos]]
                y = np.array(y)
                X += [X_arr[0][pos]]; Y += [y.mean()]; YErr += [y.std()]

            subfigMultiQT.plot(X, Y, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = "red")
            X = np.asarray(X)
            Y = np.asarray(Y)
            YErr = np.asarray(YErr)
            subfigMultiQT.fill_between(X, Y - YErr, Y + YErr, color = "red", alpha = 0.2)

            subfigMultiQT.set_xlabel(einheit_x, fontsize = labelfontsize)
            subfigMultiQT.set_ylabel('$\\chi/N$ / $J_2$', fontsize = labelfontsize)
            subfigMultiQT.set_title('Suszeptibilität pro Spin $\\chi/N$ mit bei unterschiedlichen Startvektoren', fontsize = titlefontsize)
            # plt.xticks(fontsize = 25)
            # plt.yticks(fontsize = 25)
            subfigMultiQT.tick_params(axis="both", which="major", labelsize=axisfontsize)
            # subfigMultiQT.set_xscale("log")

            legend += [Line2D([0], [0], label = "QT: " + str(max_n) + " Durchläufe", color = "blue", ls = "solid", lw = 2)]
            legend += [Line2D([0], [0], label = "QT-Mittelung", color = "red", ls = "solid", lw = 4)]

            handles, labels = plt.gca().get_legend_handles_labels()
            handles.extend(legend)

            subfigMultiQT.axhline(0, color = "grey")
            subfigMultiQT.legend(handles = handles, loc = 'best' ,frameon = False, fontsize = legendfontsize)

            for end in ends:
                subfigMultiQT.set_xlim(0.0, end)
                figMultiQT.savefig("results/" + N +  "/X_QT_N_" + N + "_J_" + J + "_0.0_" + str(end) + ".png")
                figMultiQT.savefig("results/" + N +  "/X_QT_N_" + N + "_J_" + J + "_0.0_" + str(end) + ".pdf")
                
        except:
            print("could not plot multiple runs (QT) N = %s, J = %s" %(N, J))
        plt.cla()
        plt.clf()
        plt.close(figMultiQT)
        plt.close('all')
        del figMultiQT, subfigMultiQT
        gc.collect()

import matplotlib.pyplot as plt
import numpy as np
import os
import gc
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]
N_color = [("14", "brown"), ("16", "purple")]#, ("18", "tomato")]

colors = ["red", "blue", "green", "magenta", "brown", "purple", "tomato", "cyan"]

max_n = 5

einheit_x = r'$k_B T$ / $J_2$' #'$T$ in $k_B$ / $J_2$'

titlefontsize = 40
labelfontsize = 30
legendfontsize = 30
axisfontsize = 25

start = 0.0

print("plotting specific heat (constant J1/J2, funtion of T) ...")

for end in [0.1, 0.15, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0]: # 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0
    for N_outer in range(len(N_color)):
        continue
        # if N_outer != 0: continue
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
        fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
        X_high_T = np.linspace(0.01, end, 5000)
        Y_high_T = 3/2 / 4 / X_high_T**2
        used_N = "N"
        fig1_y_max = 0; fig2_y_max = 0
        for N_inner in range(N_outer, len(N_color)):
            fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
            print("end: %f, N_outer: %s, N_inner: %s" % (end, N_color[N_outer][0], N_color[N_inner][0]))
            # ED   
            N, c = N_color[N_inner]
            file = open("results/" + N + "/data/data_specific_heat_J_const.txt", 'r')
            lines = file.readlines()
            linesJ = lines[0][len("J1/J2: "):-1]
            linesh = lines[2][len("h: "):-1]
            lbl = "N = " + N
            X = []
            Y = []
            for i in range(9,len(lines)):
                x, y = lines[i].split("\t")
                #if float(x) <= start or float(x) > end: continue
                if float(x) <= 0.0: continue
                X += [1/float(x)]
                Y += [float(y)]
            file.close()
            if max(Y) > fig1_y_max: fig1_y_max = max(Y)
            subfig1.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = lbl)
            subfig3.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = "ED")

            # QT
            file = open("results/" + N + "/data/1_data_specific_heat_J_const_QT.txt", 'r')
            lines = file.readlines()
            linesJ = lines[1][len("J1/J2: "):-1]
            linesh = lines[2][len("h: "):-1]
            lbl = "N = " + N
            X = []
            Y = []
            YErr = []
            for i in range(7,len(lines)):
                x, y, yErr = lines[i].split("\t")
                #if float(x) <= start or float(x) > end: continue
                if float(x) <= 0.0: continue
                X += [1/float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            file.close()
            subfig3.plot(X, Y, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = c, label = "QT")
            X = np.asarray(X)
            Y = np.asarray(Y)
            YErr = np.asarray(YErr)
            subfig3.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.20)

            # subfig3.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 40)
            subfig3.set_xlabel(einheit_x, fontsize = labelfontsize)
            subfig3.set_ylabel(r'$C/N$ in $J_2$', fontsize = labelfontsize)
            subfig3.set_title(r"spezifische Wärmekapazität pro Spin $C/N$ bei N = " + N, fontsize = titlefontsize)

            subfig3.axhline(0, color = "grey")
            subfig3.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

            # plt.xticks(fontsize = 25)
            # plt.yticks(fontsize = 25)
            subfig3.tick_params(axis="both", which="major", labelsize=axisfontsize)

            subfig3.set_xlim(start, end)
            # subfig3.set_ylim(0.0, 0.5)
            # subfig3.set_xscale("log")

            fig3.savefig("results/" + "C_ED_QT_" + N + "_J" + linesJ + "_h" + linesh + "_" + str(start) + "_" + str(end) + ".png")
            fig3.savefig("results/" + "C_ED_QT_" + N + "_J" + linesJ + "_h" + linesh + "_" + str(start) + "_" + str(end) + ".pdf")

            if max(Y) > fig2_y_max: fig2_y_max = max(Y)

            subfig2.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = lbl)
            subfig2.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.15)

            used_N += "_" + N

        # subfig1.set_xlabel(r'$\beta$ in $k_B$ / $J_2$', fontsize = 40)
        subfig1.set_xlabel(einheit_x, fontsize = labelfontsize)
        subfig1.set_ylabel(r'$C/N$ in $J_2$', fontsize = labelfontsize)
        subfig1.set_title(r"spezifische Wärmekapazität pro Spin $C/N$", fontsize = titlefontsize)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

        subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)
        #subfig1.tick_params(axis="both", which="major", labelsize=25)

        subfig1.set_xlim(start, end)
        # subfig1.set_ylim(0.0, fig1_y_max + 0.025)
        # subfig1.set_xscale("log")

        fig1.savefig("results/" + "C_ED_" + used_N + "_J" + linesJ + "_h" + linesh + "_" + str(start) + "_" + str(end) + ".png")
        fig1.savefig("results/" + "C_ED_" + used_N + "_J" + linesJ + "_h" + linesh + "_" + str(start) + "_" + str(end) + ".pdf")

        subfig1.plot(X_high_T, Y_high_T, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = "black", label = "high T")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        fig1.savefig("results/highT/" + "C_ED_high_T_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
        fig1.savefig("results/highT/" + "C_ED_high_T_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")


        #subfig2.set_xlabel(r'$\beta$ in $k_B$ / $J_2$', fontsize = 40)
        subfig2.set_xlabel(einheit_x, fontsize = labelfontsize)
        subfig2.set_ylabel(r'$C/N$ in $J_2$', fontsize = labelfontsize)
        subfig2.set_title(r"spezifische Wärmekapazität pro Spin $C/N$ mit einem Startvektor", fontsize = titlefontsize)

        subfig2.axhline(0, color = "grey")
        subfig2.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

        # plt.xticks(fontsize = 25)
        # plt.yticks(fontsize = 25)
        subfig2.tick_params(axis="both", which="major", labelsize=axisfontsize)

        subfig2.set_xlim(min(X), end)
        # subfig2.set_ylim(0.0, fig2_y_max + 0.025)
        # subfig2.set_xscale("log")

        fig2.savefig("results/" + "C_QT_" + used_N + "_J" + linesJ + "_h" + linesh + "_" + str(start) + "_" + str(end) + ".png")
        fig2.savefig("results/" + "C_QT_" + used_N + "_J" + linesJ + "_h" + linesh + "_" + str(start) + "_" + str(end) + ".pdf")

        subfig2.plot(X_high_T, Y_high_T, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = "black", label = "high T")
        subfig2.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
        fig2.savefig("results/highT/" + "C_QT_hight_T_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
        fig2.savefig("results/highT/" + "C_QT_hight_T_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".pdf")
        
        #plt.cla()
        plt.clf()
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)
        plt.close('all')
        del fig1, subfig1
        del fig2, subfig2
        del fig3, subfig3
        gc.collect()

    # continue

    for N_outer in range(len(N_color)):

        print("multiple runs (QT)")
        N, C = N_color[N_outer]
        for filename in os.listdir("results/" + N + "/data/excitation_energies_data/1/"):
            # continue
            if "ED" in filename: continue
            if "placeholder" in filename: continue
            figMultiQT, subfigMultiQT = plt.subplots(1,1,figsize=(16,9))
            J = filename[len("C_J"):-len("QT.txt")]
            print("end: %f, N: %s, J: %s" % (end, N_color[N_outer][0], str(J)))
            color_count = 0
            x_min = 42069
            try:
                filenameED = filename[:-len("QT.txt")] + "ED.txt"
                file = open("results/" + N + "/data/excitation_energies_data/1/" + filenameED, 'r')
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
                subfigMultiQT.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = "black", label = "ED")
                if min(X) < x_min: x_min = min(X)
            except:
                print("could not plot multiple runs (ED ref) N = %s" %N)
            try:
                for n in range(1, max_n+1):
                    file = open("results/" + N + "/data/excitation_energies_data/" + str(n) + "/" + filename, 'r')
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
                    subfigMultiQT.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = colors[color_count])
                    color_count += 1
                    if min(X) < x_min: x_min = min(X)

                subfigMultiQT.set_xlabel(einheit_x, fontsize = labelfontsize)
                subfigMultiQT.set_ylabel(r'$C/N$ in $J_2$', fontsize = labelfontsize)
                subfigMultiQT.set_title(r"spezifische Wärmekapazität pro Spin $C/N$" + "\n" + "bei unterschiedlichen Startvektoren", fontsize = titlefontsize)

                # subfigMultiQT.set_xscale("log")
                subfigMultiQT.set_xlim(x_min, end)

                # plt.xticks(fontsize = 25)
                # plt.yticks(fontsize = 25)
                subfigMultiQT.tick_params(axis="both", which="major", labelsize=axisfontsize)

                subfigMultiQT.axhline(0, color = "grey")
                # subfigMultiQT.legend(loc = 'best' ,frameon = False, fontsize = 30)
                figMultiQT.savefig("results/" + N +  "/C_QT_N_" + N + "_J_" + J + "_0.0_" + str(end) + ".png")
                figMultiQT.savefig("results/" + N +  "/C_QT_N_" + N + "_J_" + J + "_0.0_" + str(end) + ".pdf")
            except:
                print("could not plot multiple runs (QT) N = %s" %N)
            #plt.cla()
            plt.clf()
            plt.close(figMultiQT)
            plt.close('all')
            del figMultiQT, subfigMultiQT
            gc.collect()

        continue

        N, C = N_color[N_outer]
        if N_outer != 0: continue
        print("hight T for different J (ED)")
        for filename in os.listdir("results/" + N + "/data/excitation_energies_data/1/"):
            if "QT" in filename: continue
            if "placeholder" in filename: continue
            fig4, subfig4 = plt.subplots(1,1,figsize=(16,9))
            J = filename[len("C_J"):-len("ED.txt")]
            # if float(J) > 1.5: continue
            #print("%s, %s" % (N, J))
            X_high_T = np.linspace(0.01, end, 5000)
            Y_high_T = (3 * (1 + float(J)**2) / 4 / 4) / X_high_T**2
            subfig4.plot(X_high_T, Y_high_T, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = "black", label = "high T")
            used_N = "N"
            y_max = 0
            x_min = 0
            color_count = 0
            for N_inner in range(N_outer, len(N_color)):

                print("end: %f, N_outer: %s, N_inner: %s" % (end, N_color[N_outer][0], N_color[N_inner][0]))

                n, c = N_color[N_inner]

                try:
                    file = open("results/" + n + "/data/excitation_energies_data/1/" + filename, 'r')
                    lines = file.readlines()
                    X = []
                    Y = []
                    YErr = []
                    for i in range(5,len(lines)):
                        x, y = lines[i].split("\t")
                        #if float(x) < start or float(x) > end: continue
                        if float(x) <= 0.0: continue
                        if "nan" in y.lower() or "inf" in y.lower(): continue
                        X += [1/float(x)]
                        Y += [float(y)]
                    file.close()
                    subfig4.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = colors[color_count], label = "N = " + n)
                    color_count += 1
                    if min(X) < x_min: x_min = min(X)
                    if max(Y) > y_max: y_max = max(Y)
                    used_N += "_" + n
                except:
                    print("could not plot high T, N = %s and J = %s" %(n, J))

                if N_inner != len(N_color)-1: continue

                subfig4.set_xlabel(einheit_x, fontsize = labelfontsize)
                subfig4.set_ylabel(r'$C/N$ in $J_2$', fontsize = labelfontsize)
                subfig4.set_title(r"spezifische Wärmekapazität pro Spin $C/N$ bei $J_1/J_2$ = " + J, fontsize = titlefontsize)

                # subfig4.set_xscale("log")
                # subfig4.set_xlim(x_min, end)
                subfig4.set_xlim(0.0, end)
                subfig4.set_ylim(0.0, y_max + 0.025)

                # plt.xticks(fontsize = 25)
                # plt.yticks(fontsize = 25)
                subfigMultiQT.tick_params(axis="both", which="major", labelsize=axisfontsize)

                subfig4.axhline(0, color = "grey")
                subfig4.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)
                fig4.savefig("results/highT/C_J_" + J + "_" + used_N + "_0.0_" + str(end) + ".png")
                fig4.savefig("results/highT/C_J_" + J + "_" + used_N + "_0.0_" + str(end) + ".pdf")
            
            #plt.cla()
            plt.clf()
            plt.close(fig4)
            plt.close('all')
            del fig4, subfig4
            gc.collect()

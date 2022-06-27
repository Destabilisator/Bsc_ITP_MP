import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta"), ("14", "brown"), ("16", "purple"), ("18", "tomato")]

start = 0.0

print("plotting susceptibility (constant J1/J2, funtion of T) ...")

for end in [25.0, 50.0]: 
    print("from %d to %d" %(start, end))
    for N_outer in range(len(N_color)):
            fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
            fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
            used_N = "N"
            for N_inner in range(N_outer, len(N_color)):
                fig3, subfig3 = plt.subplots(1,1,figsize=(16,9))
                print("(%i,%i)" % (N_outer,N_inner))
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
                    if float(x) < start or float(x) > end: continue
                    X += [float(x)]
                    Y += [float(y)]
                file.close()
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
                    if float(x) < start or float(x) > end: continue
                    X += [float(x)]
                    Y += [float(y)]
                    YErr += [float(yErr)]
                file.close()
                subfig3.plot(X, Y, lw = 4, ls = "dashed", markersize = 0, marker = "o", color = c, label = "QT")
                X = np.asarray(X)
                Y = np.asarray(Y)
                YErr = np.asarray(YErr)
                subfig3.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.20)

                subfig3.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 40)
                subfig3.set_ylabel('$\\chi/N$ in $J_2$', fontsize = 40)
                subfig3.set_title("Suszeptibilität pro Spin $\\chi/N$ bei N = " + N, fontsize = 40)

                subfig3.axhline(0, color = "grey")
                subfig3.legend(loc = 'best' ,frameon = False, fontsize = 30)

                plt.xticks(fontsize = 25)
                plt.yticks(fontsize = 25)

                fig3.savefig("results/" + "x_ED_QT_" + N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
                plt.close(fig3)

                subfig2.plot(X, Y, lw = 4, ls = "solid", markersize = 0, marker = "o", color = c, label = lbl)
                subfig2.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.15)


                used_N += "_" + N

            subfig1.set_xlabel(r'$\beta$ in $k_B$ / $J_2$', fontsize = 40)
            subfig1.set_ylabel('$\\chi/N$ in $J_2$', fontsize = 40)
            subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$', fontsize = 40)

            subfig1.axhline(0, color = "grey")
            subfig1.legend(loc = 'best' ,frameon = False, fontsize = 30)

            # plt.xticks(fontsize = 25)
            # plt.yticks(fontsize = 25)
            subfig1.tick_params(axis="both", which="major", labelsize=25)
            #subfig1.tick_params(axis="both", which="major", labelsize=25)

            fig1.savefig("results/" + "X_ED_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.close(fig1)


            subfig2.set_xlabel(r'$\beta$ in $k_B$ / $J_2$', fontsize = 40)
            subfig2.set_ylabel('$\\chi/N$ in $J_2$', fontsize = 40)
            subfig2.set_title('Suszeptibilität pro Spin $\\chi/N$ mit einem Startvektor', fontsize = 40)

            subfig2.axhline(0, color = "grey")
            subfig2.legend(loc = 'best' ,frameon = False, fontsize = 30)

            plt.xticks(fontsize = 25)
            plt.yticks(fontsize = 25)

            fig2.savefig("results/" + "X_QT_" + used_N + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")
            plt.close(fig2)

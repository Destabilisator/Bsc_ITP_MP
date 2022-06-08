import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True

N_color = [("12", "blue")]#, ("16", "red"), ("20", "green"), ("24", "orange")] #[("6", "red"), ("8", "blue"), ("10", "green"), ("12", "orange")]#, ("14", "brown"), ("16", "purple")]

print("plotting specific heat (constant J1/J2, funtion of T) ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

used_N = "N"

for N, c in N_color:
    file = open("results/" + N + "/data/data_specific_heat_J_const.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2: "):-1]
    linesh = lines[2][len("h: "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(9,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 1, marker = "o", color = c, label = lbl)

    file = open("results/" + N + "_data_specific_heat_J_const_QT.txt", 'r')
    lines = file.readlines()
    linesJ = lines[1][len("J1/J2: "):-1]
    linesh = lines[2][len("h: "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    YErr = []
    for i in range(7,len(lines)):
        x, y, yErr = lines[i].split("\t")
        # if float(x) > 0:
        X += [float(x)]
        Y += [float(y)]
        YErr += [float(yErr)]
    subfig1.plot(X, Y, lw = 1, ls = "dashed", markersize = 1, marker = "o", color = c, label = "QT", alpha = 0.2)
    X = np.asarray(X)
    Y = np.asarray(Y)
    YErr = np.asarray(YErr)
    subfig1.fill_between(X, Y - YErr, Y + YErr, color = c, alpha = 0.05)

    used_N += "_" + N

subfig1.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 18)
# subfig1.set_xlabel(r'$T$ in $k_B$ / $J_2$', fontsize = 18)
subfig1.set_ylabel(r'spezifische Wärmekapazität pro Spin $C/N$ in $J_2$', fontsize = 18)
subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)

subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "C_" used_N + "_J" + linesJ + "_h" + linesh + ".png")
#plt.savefig("results/" + "specific_heat_J_const.png")
#plt.show

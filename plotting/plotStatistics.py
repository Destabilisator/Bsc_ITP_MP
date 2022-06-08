import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "magenta")]#[("16", "red"), ("18", "blue"), ("20", "green"), ("22", "magenta"), ("24", "brown"), ("26", "purple"), ("28", "tomato")]
n_color = [("1", "red"), ("2", "blue"), ("3", "green"), ("4", "tomato")]

def plot_n_for_each_N(start: float, end: float):
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
            if float(x) < start or float(x) > end: continue
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
                if float(x) < start or float(x) > end: continue
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
        plt.savefig("results/QT_stats/C_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")

def plot_n_for_each_N_sigma(start: float, end: float): 
    print("plotting dependance of n for fixed N (plotting only sigma) ...")
    for N, NC in N_color:
        print("N = " + N)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
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
                if float(x) < start or float(x) > end: continue
                X += [float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            subfig1.plot(X, YErr, lw = 1, ls = "solid", markersize = 0, marker = "o", color = nc, label = "QT: n = " + n)
        # saving
        subfig1.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
        subfig1.set_ylabel(r'$\sigma$ in $J_2$', fontsize = 25)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        subfig1.set_title(r"Standardabweichung $\sigma$ von $C/N$ bei der QT mit N = " + N + r" und $J_1$ / $J_2$ = " + linesJ, fontsize = 25)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)
        filename = "N_" + N + "_n"
        for n, nc in n_color:
            filename += "_" + n
        plt.savefig("results/QT_stats/C_sigma_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")

def plot_N_for_each_n(start: float, end: float):
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
                if float(x) < start or float(x) > end: continue
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
                if float(x) < start or float(x) > end: continue
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
        plt.savefig("results/QT_stats/C_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")

def plot_N_for_each_n_sigma(start: float, end: float):
    print("plotting dependance of N for fixed n (plotting only sigma) ...")
    for n, nc in n_color:
        print("n = " + n)
        fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
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
                if float(x) < start or float(x) > end: continue
                X += [float(x)]
                Y += [float(y)]
                YErr += [float(yErr)]
            subfig1.plot(X, YErr, lw = 1, ls = "solid", markersize = 0, marker = "o", color = NC, label = "QT: N = " + N, alpha = 1.0)
        # saving
        subfig1.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
        subfig1.set_ylabel(r'$\sigma$ in $J_2$', fontsize = 25)
        #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
        subfig1.set_title(r"Standardabweichung $\sigma$ von $C/N$ mit $J_1$ / $J_2$ = " + linesJ + " nach " + n + " Startvektoren", fontsize = 25)
        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)
        filename = "n_" + n + "_N"
        for N, NC in N_color:
            filename += "_" + N
        plt.savefig("results/QT_stats/C_sigma_" + filename + "_J" + linesJ + "_" + str(start) + "_" + str(end) + ".png")

def plot_delta_ED(start: float, end: float):
    print("plotting \Delta E ED - QT for fixed N ...")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    used_N = "N"
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
            if float(x) < start or float(x) > end: continue
            X_QT += [float(x)]
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
            if float(x) < start or float(x) > end: continue
            X_ED += [float(x)]
            Y_ED += [float(y)]

        X = []
        Y = []

        # print(len(X_QT))
        # print(len(Y_QT))
        # print(len(X_ED))
        # print(len(Y_ED))

        for i in range(0,len(X_QT)):
            X += [X_QT[i]] # [(X_QT[i] + X_ED[i])/2]
            Y += [Y_QT[i] - Y_ED[i]]

        subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 0, marker = "o", color = NC, label = "QT-ED: N = " + N, alpha = 1.0)

        used_N += "_" + N

    # saving
    subfig1.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
    subfig1.set_ylabel(r'Abweichung in $J_2$', fontsize = 25)
    #subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 18)
    subfig1.set_title(r"Abweichung von $C/N$ zwischen ED und QT mit $J_1$ / $J_2$ = " + linesJQT, fontsize = 25)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)
    plt.savefig("results/QT_stats/C_delta_ED_QT_" + used_N + "_J" + linesJQT + "_" + str(start) + "_" + str(end) + ".png")

if __name__ == "__main__":
    print("plotting specific heat (constant J1/J2, funtion of T):")
    start = float(sys.argv[1])
    end = float(sys.argv[2])
    plot_n_for_each_N(start, end)
    print()
    plot_n_for_each_N_sigma(start, end)
    print()
    plot_N_for_each_n(start, end)
    print()
    plot_N_for_each_n_sigma(start, end)
    print()
    plot_delta_ED(start, end)
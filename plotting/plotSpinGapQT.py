import matplotlib.pyplot as plt
import os
import sys
import numpy as np
plt.rcParams['text.usetex'] = True

N_color = []
N_color_LOW = [("10", "green"), ("12", "magenta"), ("14", "brown")]#, ("16", "purple"), ("18", "tomato")] # ("6", "red"), ("8", "blue"), 
N_color_HIGH = [("20", "red"), ("22", "blue"), ("24", "green"), ("26", "magenta")]#, ("28", "brown"), ("30", "purple"), ("32", "tomato")]
no_ED = False

peak_offset = 2000 #1500 # 1000
fit_samples = 10
search_start_percent = 4/5
search_end_percent = 1/5

def sort_data(X, Y):
    length = len(X)
    for i in range(length-1):
        for j in range(0, length-i-1):
            x1 = X[j]; x2 = X[j+1]
            if x1 > x2:
                X[j], X[j+1] = X[j+1], X[j]
                Y[j], Y[j+1] = Y[j+1], Y[j]
    return X, Y

def expFunc(x: float, A: float, k: float, x_0: float, y_0: float) -> float:
    return A * np.exp(k * (x + x_0)) + y_0

def get_spin_gap(N, J, filename) -> (float, float, float, float):
    file = open("results/" + N + "/data/spin_gap_data/" + filename, 'r')
    ED_QT = filename[len(filename)-6:-4]
    lines = file.readlines()
    fig2, subfig2 = plt.subplots(1,1,figsize=(16,9))
    X = []; Y = []
    for i in range(5,len(lines)):
        arr = lines[i].split("\t")
        x = arr[0]
        y = arr[1]
        X += [float(x)] # beta -> T
        Y += [float(y)]

    X.reverse(); Y.reverse()

    X_fit = X[0:int(len(X)*search_start_percent)]; Y_fit = Y[0:int(len(X)*search_start_percent)]
    X_fit = np.array(X_fit); Y_fit = np.array(Y_fit)
    fitting_results = np.polyfit(X_fit, Y_fit, 1, full = True)
    m, b = fitting_results[0]
    SSE = fitting_results[1][0]
    diff = Y_fit - Y_fit.mean()
    square_diff = diff ** 2
    SST = square_diff.sum()
    R2 = 1 - SSE/SST
    X_fit_range = X_fit; Y_fit_range = Y_fit
    X_fit = np.array(X_fit); Y_fit = np.array(Y_fit)

    for p in range(1, fit_samples+1):
        percent = search_start_percent + (search_end_percent - search_start_percent) * float(p) / float(fit_samples)
        X_fit = []; Y_fit = []
        for i in range(0, int(len(X)*percent)):
            X_fit += [X[i]]; Y_fit += [np.log(Y[i])] # log = ln [np.log(Y[i])]
        X_fit = np.array(X_fit); Y_fit = np.array(Y_fit)
        fitting_results = np.polyfit(X_fit, Y_fit, 1, full = True)
        m_new, b_new = fitting_results[0]
        SSE = fitting_results[1][0]
        diff = Y_fit - Y_fit.mean()
        square_diff = diff ** 2
        SST = square_diff.sum()
        R2_new = 1 - SSE/SST
        if R2_new > R2:
            R2 = R2_new
            m = m_new
            b = b_new
            X_fit_range = X_fit; Y_fit_range = Y_fit

    Y_fitted_range = np.exp(Y_fit_range)
    subfig2.plot(X_fit_range, Y_fitted_range, lw = 1, ls = "solid", markersize = 5, marker = "o", color = "green", label = "range")

    #print("%f, %f, %f" % (R2, m, b))

    # X_fit = []; Y_fit = []
    # for i in range(0, len(X)-peak_offset):
    #     if Y[i] > Y[i+peak_offset]: break
    #     X_fit += [X[i]]; Y_fit += [np.log(Y[i])] # log = ln [np.log(Y[i])]
    # if (len(X_fit) == 0):
    #     for i in range(0, len(X) - int(len(X)/3)):
    #         X_fit += [X[i]]; Y_fit += [np.log(Y[i])]
    # X_fit = np.array(X_fit); Y_fit = np.array(Y_fit)
    # # fitting
    # fitting_results = np.polyfit(X_fit, Y_fit, 1, full = True)
    # # evaluate quality of fit
    # m, b = fitting_results[0]
    # SSE = fitting_results[1][0]
    # diff = Y_fit - Y_fit.mean()
    # square_diff = diff ** 2
    # SST = square_diff.sum()
    # R2 = 1 - SSE/SST
    # print(R2)

    subfig2.plot(X, Y, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "blue", label = "QT data")
    Y_fitted = [np.exp(m * x + b) for x in X_fit_range]
    subfig2.plot(X_fit_range, Y_fitted, lw = 1, ls = "solid", markersize = 1, marker = "o", color = "red", label = "exp fit, R = " + str(R2))
    #subfig2.set_xlabel(r'T in $k_B$ / $J_2$', fontsize = 25)
    subfig2.set_xlabel(r'$\beta$ in $J_2$ / $k_B$', fontsize = 25)
    subfig2.set_ylabel('$\\chi/N$ in $J_2$', fontsize = 25)
    subfig2.set_title("$\\chi/N$ für N = " + N + r" mit $J_1$ / $J_2$ = " + J, fontsize = 25)

    # subfig2.axhline(0, color = "grey")
    subfig2.legend(loc = 'best' ,frameon = False, fontsize = 20)

    plt.axvline(x=X_fit_range[len(X_fit_range)-1], color='black', linestyle='--')
    #subfig2.set_xscale("log")
    subfig2.set_yscale("log")
    plt.savefig("results/" + N + "/data/spin_gap_data/X_J" + J + "_" + ED_QT + ".png")
    plt.close(fig2)
    return np.exp(b), abs(m), 0.0, 0.0


if __name__ == "__main__":

    if len(sys.argv) > 1:
        regime = sys.argv[1]
        if regime == "low": N_color = N_color_LOW#; print("low regime")
        elif regime == "high": N_color = N_color_HIGH; no_ED = True#; print("high regime")
        else: N_color = N_color_LOW; print("default low (wrong args)")
    else: N_color = N_color_LOW; print("default low (no args)")

    print("plotting spin gap ...")
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

    used_N = "N"

    for N, c in N_color:
        print("N = " + N + ":") 
        # QT results
        print("QT (exp fit)...")
        #lbl = "QT (exp fit): N = " + N
        lbl = "N = " + N
        X = []; Y = []
        for filename in os.listdir("results/" + N + "/data/spin_gap_data/"):
            if filename[len(filename)-6:] != "QT.txt": continue
            J = filename[len("X_J"):-len("QT.txt")]
            # print(J)
            A, k, x_0, y_0 = get_spin_gap(N, J, filename)
            # print(str(A) + " " + str(k)  + " " + str(x_0)  + " " + str(y_0))
            X += [float(J)]
            Y += [float(k)]
        X, Y = sort_data(X, Y)
        subfig1.plot(X, Y, lw = 1, ls = "dashed", markersize = 0, marker = "o", color = c, label = lbl)

        if not no_ED:
            # ED results exp fit
            print("ED (exp fit)...")
            lbl = "ED (exp fit): N = " + N
            X = []; Y = []
            for filename in os.listdir("results/" + N + "/data/spin_gap_data/"):
                if filename[len(filename)-6:] != "ED.txt": continue
                J = filename[len("X_J"):-len("ED.txt")]
                # print(J)
                A, k, x_0, y_0 = get_spin_gap(N, J, filename)
                # print(str(A) + " " + str(k)  + " " + str(x_0)  + " " + str(y_0))
                X += [float(J)]
                Y += [float(k)]
            X, Y = sort_data(X, Y)
            subfig1.plot(X, Y, lw = 0, ls = "dotted", markersize = 2, marker = "o", color = c, alpha = 0.5)#, label = lbl, alpha = 0.5)
            # ED results
            print("ED (dispersion)...")
            lbl = "ED : N = " + N
            X = []; Y = []
            file = open("results/" + N + "/data/data_spin_gap.txt", 'r') # _data_spin_gap / _data_spin_gap_with_index
            lines = file.readlines()
            lbl = "ED: N = " + N
            X = []; Y = []
            for i in range(7,len(lines)):
                arr = lines[i].split("\t")
                x = arr[0]
                y = arr[1]
                X += [float(x)]
                Y += [float(y)]
            subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 0, marker = "o", color = c, alpha = 0.4) #  label = lbl,
        
        print()

        used_N += "_" + N

        subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
        subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 25)
        subfig1.set_title(r'Spingap Energien $\Delta E_{gap}$', fontsize = 25)
        # subfig1.set_title(r'Spingap Energies $\Delta E_{gap}$ für $\Delta = 1$', fontsize = 18)

        subfig1.axhline(0, color = "grey")
        subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

        plt.savefig("results/" + "spin_gap_with_QT_" + used_N + ".png")
        #plt.show()

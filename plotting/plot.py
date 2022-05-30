import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
plt.rcParams['text.usetex'] = True
import sys
import time

def plot_delta_E(N):
    print("plotting delta E ...")
    file = open("results/" + N + "_data_delta_E.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + N + "\n" +"\\# Datenpunkte: " + lines[3][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(7,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 18)
    subfig1.set_title(r'Anregungsenergieren $\Delta E$', fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_delta_E.png")

def plot_specific_heat_T_const(N):
    print("plotting specific heat (constant T, funtion of J1/J2) ...")
    file = open("results/" + N + "_data_specific_heat_T_const.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("T = "):-1]
    lbl = r"$T$ = " + linesBeta + "\n" + "N = " + N + "\n" + "\\# Datenpunkte: " + lines[4][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel(r'spezifische Wärmekapazität pro Spin $C/N$ in $J_2$', fontsize = 18)
    subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $T$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_specific_heat_T_const.png")

def plot_specific_heat_J_const(N):
    print("plotting specific heat (constant J1/J2, funtion of T) ...")
    file = open("results/" + N + "_data_specific_heat_J_const.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$T$ $k_B$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel(r'spezifische Wärmekapazität pro Spin $C/N$ in $J_2$', fontsize = 18)
    subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", $k_B$ = 1", fontsize = 18)

    file = open("results/" + N + "_data_specific_heat_J_const_QT.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(3,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'red', label = "QT")

    file = open("results/" + N + "_data_specific_heat_J_const_QT_multi.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'green', label = "QT")

    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_specific_heat_J_const.png")

def plot_susceptibility_T_const(N):
    print("plotting suszeptibility (constant T, funtion of J1/J2) ...")
    file = open("results/" + N + "_data_susceptibility_T_const.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("T = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel('Suszeptibilität pro Spin $\\chi/N$ in $J_2$', fontsize = 18)
    subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$ mit $T$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)
    #subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_susceptibility_T_const.png")

def plot_susceptibility_J_const(N):
    print("plotting suszeptibility (constant J1/J2, funtion of T) ...")
    file = open("results/" + N + "_data_susceptibility_J_const.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$T$ $k_B$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel('Suszeptibilität pro Spin $\\chi/N$ in $J_2$', fontsize = 18)
    subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$ für ' + r"$J_1$ / $J_2$ = " + linesJ + r", $k_B$ = 1", fontsize = 18)
    

    # file = open("results/" + N + "_momentum_susceptibility_J_const.txt", 'r')
    # lines = file.readlines()
    # linesJ = lines[0][len("J1/J2 = "):-1]
    # lbl = "N = " + N
    # X = []
    # Y = []
    # for i in range(8,len(lines)):
    #     x, y = lines[i].split("\t")
    #     X += [float(x)]
    #     Y += [float(y)]
    # subfig1.plot(X, Y, lw = 0.5, ls = "solid", markersize = 1, marker = "o", color = 'red', label = "MS")

    # file = open("results/" + N + "_spininversion_susceptibility_J_const.txt", 'r')
    # lines = file.readlines()
    # linesJ = lines[0][len("J1/J2 = "):-1]
    # lbl = "N = " + N
    # X = []
    # Y = []
    # for i in range(8,len(lines)):
    #     x, y = lines[i].split("\t")
    #     X += [float(x)]
    #     Y += [float(y)]
    # subfig1.plot(X, Y, lw = 0.5, ls = "solid", markersize = 1, marker = "o", color = 'green', label = "SI")

    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_susceptibility_J_const.png")

def plot_k_dispersion_J_const(N):
    print("plotting energy dispersion (constant J1/J2) ...")
    file = open("results/" + N + "_momentum_energy_dispersion_J_const.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = 500, color = 'blue', label = lbl)
    plt.ticklabel_format(style='plain', axis='x', useOffset=False)
    plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
    subfig1.set_xlabel(r'Impulsquantenzahl $k$', fontsize = 18)
    subfig1.set_ylabel('Energie E in $J_2$', fontsize = 18)
    subfig1.set_title('Energiedispersion für ' + r"$J_1$ / $J_2$ = " + linesJ, fontsize = 18)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_energy_dispersion_J_const.png")

def plot_spin_gap(N):
    print("plotting spin gap ...")
    file = open("results/" + N + "_data_spin_gap.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + N + "\n" +"\\# Datenpunkte: " + lines[3][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(7,len(lines)):
        x, y = lines[i].split("\t")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 18)
    subfig1.set_title(r'Spingap Energies $\Delta E_{gap}$ für $\Delta = 1$', fontsize = 18)
    subfig1.axhline(0, color = "grey")
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/" + N + "_spin_gap.png")

if __name__ == "__main__":

    start_time = time.time()

    N = sys.argv[1]
    
    # plot_delta_E(N)
    # plot_specific_heat_T_const(N)
    plot_specific_heat_J_const(N)
    # if sys.argv[3] != "noX":
    #     plot_susceptibility_T_const(N)
    # plot_susceptibility_J_const(N)
    # plot_k_dispersion_J_const(N)
    # plot_spin_gap(N)

    end_time = time.time()

    print("plotting done. This took %.2f seconds" % (end_time - start_time) )

    if (sys.argv[2] == "show"):
        plt.show()
        
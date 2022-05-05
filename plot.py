import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import sys

def plot_delta_E(N):
    print("plotting delta E, N = " + N + " ...")
    file = open("results/" + N + "_data_delta_E.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + N + "\n" +"\\# Datenpunkte: " + lines[3][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(7,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 18)
    subfig1.set_title(r'Anregungsenergieren $\Delta E$', fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_data_delta_E.png")

def plot_specific_heat(N):
    print("plotting specific heat (constant T, funtion of J1/J2), N = " + N + " ...")
    file = open("results/" + N + "_data_specific_heat.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("T = "):-1]
    lbl = r"$T$ = " + linesBeta + "\n" + "N = " + N + "\n" + "\\# Datenpunkte: " + lines[4][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
    subfig1.set_ylabel(r'spezifische Wärmekapazität pro Spin $C/N$ in $J_2$', fontsize = 18)
    subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $T$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_data_specific_heat.png")

def plot_momentum_specific_heat(N):
    print("plotting specific heat (constant J1/J2, funtion of T), N = " + N + " ...")
    file = open("results/" + N + "_momentum_specific_heat.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = r"$J_1$ / $J_2$ = " + linesJ + "\n" + "N = " + N + "\n" + "\\# Datenpunkte: " + lines[4][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$T$ in $J_2$/$k_B$', fontsize = 18)
    subfig1.set_ylabel(r'spezifische Wärmekapazität pro Spin $C/N$ in $J_2$', fontsize = 18)
    subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $J_1$ / $J_2$ = ' + linesJ + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_momentum_specific_heat.png")

def plot_naiv_susceptibility(N):
    print("plotting suszeptibility (constant J1/J2, funtion of T), N = " + N + " ...")
    file = open("results/" + N + "_naiv_susceptibility.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = r"$J_1$ / $J_2$ = " + linesJ + "\n" + "N = " + N + "\n" + "\\# Datenpunkte: " + lines[4][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$T$ in $J_2$/$k_B$', fontsize = 18)
    subfig1.set_ylabel('Suszeptibilität pro Spin $\\chi/N$ in $J_2$', fontsize = 18)
    subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$ für ' + r"$J_1$ / $J_2$ = " + linesJ + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_naiv_susceptibility.png")

def plot_magnetization_susceptibility(N):
    print("plotting suszeptibility (constant J1/J2, funtion of T), N = " + N + " ...")
    file = open("results/" + N + "_magnetization_susceptibility.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2 = "):-1]
    lbl = r"$J_1$ / $J_2$ = " + linesJ + "\n" + "N = " + N + "\n" + "\\# Datenpunkte: " + lines[4][len("datapoints: "):-1]
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'$T$ in $J_2$/$k_B$', fontsize = 18)
    subfig1.set_ylabel('Suszeptibilität pro Spin $\\chi/N$ in $J_2$', fontsize = 18)
    subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$ für ' + r"$J_1$ / $J_2$ = " + linesJ + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_magnetization_susceptibility.png")



if __name__ == "__main__":

    N = sys.argv[1]
    
    plot_delta_E(N)
    plot_specific_heat(N)
    plot_momentum_specific_heat(N)
    #plot_naiv_susceptibility(N)
    plot_magnetization_susceptibility(N)

    plt.show()
        
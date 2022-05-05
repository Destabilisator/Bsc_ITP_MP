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
    print("plotting specific heat (constant beta, funtion of J1/J2), N = " + N + " ...")
    file = open("results/" + N + "_data_specific_heat.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("beta = "):-1]
    lbl = r"$\beta$ = " + linesBeta + "\n" + "N = " + N + "\n" + "\\# Datenpunkte: " + lines[4][len("datapoints: "):-1]
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
    subfig1.set_ylabel(r'spezifische Wärmekapazität $C$  in $J_2$', fontsize = 18)
    subfig1.set_title(r'spezifische Wärmekapazitäten $C$ mit $\beta$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_data_specific_heat.png")

def plot_momentum_specific_heat(N):
    print("plotting specific heat (constant J1/J2, funtion of beta), N = " + N + " ...")
    file = open("results/" + N + "_momentum_specific_heat.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("beta = "):-1]
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
    subfig1.set_xlabel(r'$\beta$', fontsize = 18)
    subfig1.set_ylabel(r'spezifische Wärmekapazität $C$  in $J_2$', fontsize = 18)
    subfig1.set_title(r'spezifische Wärmekapazitäten $C$ mit $J_1$ / $J_2$ = ' + linesJ + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_momentum_specific_heat.png")

def plot_momentum_magnetization(N):
    print("plotting magnetization (constant J1/J2, funtion of beta), N = " + N + " ...")
    file = open("results/" + N + "_naiv_magnetization.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("beta = "):-1]
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
    subfig1.set_xlabel(r'$\beta$', fontsize = 18)
    subfig1.set_ylabel('Magnetisierung $\\chi$  in $J_2$', fontsize = 18)
    subfig1.set_title('Magnetisierung $\\chi$ für N = ' + N + r" und $J_1$ / $J_2$ = " + linesJ + r", $k_B$ = 1", fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/" + N + "_naiv_magnetization.png")



if __name__ == "__main__":

    N = sys.argv[1]
    
    plot_delta_E(N)
    plot_specific_heat(N)
    plot_momentum_specific_heat(N)
    plot_momentum_magnetization(N)

    #plt.show()
        
import matplotlib.pyplot as plt

def plot_delta_E():
    file = open("results/data_delta_E.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + lines[0][len("N: "):-1] + "\n" +"# Datenpunkte: " + lines[3][len("datapoints: "):-1]
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
    subfig1.set_title(r'Anregungsenergieren $\Delta E$ für N = ' + lines[0][len("N: "):-1], fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    subfig1.axhline(0, color = "grey")
    plt.savefig("results/data_delta_E.png")

def plot_specific_heat():
    file = open("results/data_specific_heat.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("beta = "):-1]
    linesN = lines[1][len("N: "):-1]
    lbl = r"$\beta$ = " + linesBeta + "\n" + "N = " + linesN + "\n" + "# Datenpunkte: " + lines[4][len("datapoints: "):-1]
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
    subfig1.set_title(r'spezifische Wärmekapazitäten $C$ für N = ' + linesN + r" und $\beta$ = " + linesBeta, fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/data_specific_heat.png")

def plot_momentum_specific_heat():
    file = open("results/momentum_specific_heat.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("beta = "):-1]
    linesN = lines[1][len("N: "):-1]
    lbl = r"$J_1$ / $J_2$ = " + linesJ + "\n" + "N = " + linesN + "\n" + "# Datenpunkte: " + lines[4][len("datapoints: "):-1]
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
    subfig1.set_title(r'spezifische Wärmekapazitäten $C$ für N = ' + linesN + r" und $J_1$ / $J_2$ = " + linesJ, fontsize = 18)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    #subfig1.axhline(0, color = "grey")
    plt.savefig("results/momentum_specific_heat.png")

if __name__ == "__main__":
    print("plotting delta E ...")
    plot_delta_E()
    print("plotting specific heat (constant beta, funtion of J1/J2) ...")
    plot_specific_heat()
    print("plotting specific heat (constant J1/J2, funtion of beta) ...")
    plot_momentum_specific_heat()
    plt.show()
    
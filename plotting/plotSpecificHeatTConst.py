import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "orange")]#, ("14", "brown"), ("16", "purple")]

print("plotting specific heat (constant T, funtion of J1/J2) ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

used_N = "N"

for N, c in N_color:
    file = open("results/" + N + "/data/data_specific_heat_T_const.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("beta: "):-1]
    linesh = lines[2][len("h: "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(9,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = c, label = lbl)

    used_N += "_" + N

subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
subfig1.set_ylabel(r'$C/N$ in $J_2$', fontsize = 25)
subfig1.set_title(r'$C/N$ mit $\beta$ = ' + linesBeta + r", h = " + linesh + r" und $k_B$ = 1", fontsize = 25)

#subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

plt.savefig("results/" + "C_" + used_N + "_beta" + linesBeta + "_h" + linesh + ".png")
#plt.show()

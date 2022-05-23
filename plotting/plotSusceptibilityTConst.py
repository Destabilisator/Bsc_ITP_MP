import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "orange")]#, ("14", "brown"), ("16", "purple")]

print("plotting suszeptibility (constant T, funtion of J1/J2) ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

for N, c in N_color:
    file = open("results/" + N + "_data_susceptibility_T_const.txt", 'r')
    lines = file.readlines()
    linesBeta = lines[0][len("T = "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = c, label = lbl)

subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
subfig1.set_ylabel('Suszeptibilität pro Spin $\\chi/N$ in $J_2$', fontsize = 18)
subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$ mit $T$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)

#subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "susceptibility_T_const.png")
#plt.show()

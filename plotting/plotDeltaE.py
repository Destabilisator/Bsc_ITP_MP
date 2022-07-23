import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "orange"), ("14", "brown"), ("16", "purple")]

print("plotting Delta E ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

for N, c in N_color:
    file = open("results/" + N + "/data/data_delta_E.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + N
    linesh = lines[1][len("h: "):-1]
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 3, ls = "solid", markersize = 5, marker = "o", color = c, label = lbl)

subfig1.tick_params(axis="both", which="major", labelsize=25)

subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 25)
subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 25)
subfig1.set_title(r'Anregungsenergieren $\Delta E$', fontsize = 25)

subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

plt.savefig("results/" + "delta_E.png")
#plt.show()

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "orange")]#, ("14", "brown"), ("16", "purple"), ("20", "cyan")] #, ("18", "black")

print("plotting spin gap ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

for N, c in N_color:
    file = open("results/" + N + "_data_spin_gap_with_index.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(7,len(lines)):
        arr = lines[i].split("\t")
        x = arr[0]
        y = arr[1]
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = c, label = lbl)

subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
subfig1.set_ylabel(r'$\Delta E_{gap}$  in $J_2$', fontsize = 18)
subfig1.set_title(r'Spingap Energies $\Delta E_{gap}$ für $\Delta = 1$', fontsize = 18)

subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "spin_gap.png")
#plt.show()

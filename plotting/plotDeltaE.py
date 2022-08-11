import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

titlefontsize = 39
labelfontsize = 35
legendfontsize = 35
axisfontsize = 30

N_color = [("10", "red"), ("12", "blue"), ("14", "green"), ("16", "magenta")]
N_color = [("8", "purple"), ("10", "red"), ("12", "blue"), ("14", "green"), ("16", "magenta")]
# N_color = [("6", "brown"), ("8", "purple"), ("10", "red"), ("12", "blue"), ("14", "green"), ("16", "magenta")]

print("plotting Delta E ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
used_N = "_N"
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
    used_N += "_" + N

subfig1.tick_params(axis="both", which="major", labelsize=axisfontsize)

subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = labelfontsize)
subfig1.set_ylabel(r'$\Delta E$ / $J_2$ = $(E_1 - E_0)$ / $J_2$', fontsize = labelfontsize)
subfig1.set_title(r'Anregungsenergieren $\Delta E$', fontsize = titlefontsize)

subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = legendfontsize)

plt.savefig("results/" + "delta_E" + used_N + "_.png")
#plt.show()

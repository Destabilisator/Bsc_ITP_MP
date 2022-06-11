import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
plt.rcParams['text.usetex'] = True

N_color_sz = [("6", "red", 5000), ("8", "blue", 2000), ("10", "green", 750), ("12", "orange", 100)]#, ("14", "brown", 50), ("16", "purple", 25)]

print("plotting energy dispersion (constant J1/J2) ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))

used_N = "N"

for N, c, sz in N_color_sz:
    file = open("results/" + N + "/data/momentum_energy_dispersion_J_const.txt", 'r')
    lines = file.readlines()
    linesJ = lines[0][len("J1/J2: "):-1]
    linesh = lines[2][len("h: "):-1]
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(8,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = sz, color = c, label = lbl)
    
    used_N += "_" + N

plt.ticklabel_format(style='plain', axis='x', useOffset=False)
plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
subfig1.set_xlabel(r'Impulsquantenzahl $k$', fontsize = 25)
subfig1.set_ylabel('Energie E in $J_2$', fontsize = 25)
subfig1.set_title(r"Energiedispersion bei $J_1$ / $J_2$ = " + linesJ + " und h = " + linesh + r" mit $k_B$ = 1", fontsize = 25)

subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 20)

plt.savefig("results/" + "energy_dispersion_" + used_N + "_J" linesJ ".png")
#plt.show()

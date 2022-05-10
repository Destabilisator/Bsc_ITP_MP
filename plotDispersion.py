import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
plt.rcParams['text.usetex'] = True

N = "6"

print("plotting energy dispersion (constant J1/J2) ...")
file = open("results/" + N + "_momentum_energy_dispersion_J_const.txt", 'r')
lines = file.readlines()
linesJ = lines[0][len("J1/J2 = "):-1]
lbl = "N = " + N
X = []
Y = []
for i in range(8,len(lines)):
    x, y = lines[i].split("\t")
    #print(x + " " + y + "\r")
    X += [float(x)]
    Y += [float(y)]
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = 5000, color = 'red', label = lbl)
plt.ticklabel_format(style='plain', axis='x', useOffset=False)
plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
subfig1.set_xlabel(r'Impulsquantenzahl $k$', fontsize = 18)
subfig1.set_ylabel('Energie E in $J_2$', fontsize = 18)
subfig1.set_title('Energiedispersion für ' + r"$J_1$ / $J_2$ = " + linesJ, fontsize = 18)


N = "8"

file = open("results/" + N + "_momentum_energy_dispersion_J_const.txt", 'r')
lines = file.readlines()
linesJ = lines[0][len("J1/J2 = "):-1]
lbl = "N = " + N
X = []
Y = []
for i in range(8,len(lines)):
    x, y = lines[i].split("\t")
    #print(x + " " + y + "\r")
    X += [float(x)]
    Y += [float(y)]
subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = 2000, color = 'blue', label = lbl)

N = "10"

file = open("results/" + N + "_momentum_energy_dispersion_J_const.txt", 'r')
lines = file.readlines()
linesJ = lines[0][len("J1/J2 = "):-1]
lbl = "N = " + N
X = []
Y = []
for i in range(8,len(lines)):
    x, y = lines[i].split("\t")
    #print(x + " " + y + "\r")
    X += [float(x)]
    Y += [float(y)]
subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = 750, color = 'green', label = lbl)

N = "12"

file = open("results/" + N + "_momentum_energy_dispersion_J_const.txt", 'r')
lines = file.readlines()
linesJ = lines[0][len("J1/J2 = "):-1]
lbl = "N = " + N
X = []
Y = []
for i in range(8,len(lines)):
    x, y = lines[i].split("\t")
    #print(x + " " + y + "\r")
    X += [float(x)]
    Y += [float(y)]
subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = 100, color = 'orange', label = lbl)

# N = "14"

# file = open("results/" + N + "_momentum_energy_dispersion_J_const.txt", 'r')
# lines = file.readlines()
# linesJ = lines[0][len("J1/J2 = "):-1]
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(8,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# subfig1.scatter(X, Y, linewidth = 1, marker = "_", s = 50, color = 'brown', label = lbl)

subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "energy_dispersion_J_const.png")
#plt.show()

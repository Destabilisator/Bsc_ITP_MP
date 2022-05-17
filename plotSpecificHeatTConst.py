import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

N = "6"

print("plotting specific heat (constant T, funtion of J1/J2) ...")
file = open("results/" + N + "_data_specific_heat.txt", 'r')
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
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'red', label = lbl)
subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
subfig1.set_ylabel(r'spezifische Wärmekapazität pro Spin $C/N$ in $J_2$', fontsize = 18)
subfig1.set_title(r'spezifische Wärmekapazität pro Spin $C/N$ mit $T$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)

N = "8"

file = open("results/" + N + "_data_specific_heat.txt", 'r')
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
subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)

N = "10"

file = open("results/" + N + "_data_specific_heat.txt", 'r')
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
subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'green', label = lbl)

N = "12"

file = open("results/" + N + "_data_specific_heat.txt", 'r')
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
subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'orange', label = lbl)

N = "14"

file = open("results/" + N + "_data_specific_heat.txt", 'r')
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
subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'brown', label = lbl)

# N = "16"

# file = open("results/" + N + "_data_specific_heat.txt", 'r')
# lines = file.readlines()
# linesBeta = lines[0][len("T = "):-1]
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(8,len(lines)):
#    x, y = lines[i].split("\t")
#    #print(x + " " + y + "\r")
#    X += [float(x)]
#    Y += [float(y)]
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'purple', label = lbl)

subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "specific_heat_T_const.png")
#plt.show()

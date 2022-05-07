import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

N = "6"

print("plotting suszeptibility (constant J1/J2, funtion of T) ...")
file = open("results/" + N + "_data_magnetic_susceptibility_T_const.txt", 'r')
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
subfig1.set_ylabel('Suszeptibilität pro Spin $\\chi/N$ in $J_2$', fontsize = 18)
subfig1.set_title('Suszeptibilität pro Spin $\\chi/N$ mit $T$ = ' + linesBeta + r", $k_B$ = 1", fontsize = 18)
subfig1.axhline(0, color = "grey")

N = "8"

file = open("results/" + N + "_data_magnetic_susceptibility_T_const.txt", 'r')
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

file = open("results/" + N + "_data_magnetic_susceptibility_T_const.txt", 'r')
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

file = open("results/" + N + "_data_magnetic_susceptibility_T_const.txt", 'r')
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

#N = "14"

#print("plotting suszeptibility (constant J1/J2, funtion of T) ...")
#file = open("results/" + N + "_data_magnetic_susceptibility_T_const.txt", 'r')
#lines = file.readlines()
#linesBeta = lines[0][len("T = "):-1]
#lbl = "N = " + N
#X = []
#Y = []
#for i in range(8,len(lines)):
#    x, y = lines[i].split("\t")
#    #print(x + " " + y + "\r")
#    X += [float(x)]
#    Y += [float(y)]
#subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'brown', label = lbl)


subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "susceptibility_T_const.png")
#plt.show()

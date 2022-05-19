import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


N_color = [("6", "red"), ("8", "blue"), ("10", "green"), ("12", "orange")]#, ("14", "brown"), ("16", "purple")]

print("plotting Delta E ...")
fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 18)
subfig1.set_title(r'Anregungsenergieren $\Delta E$', fontsize = 18)

for N, c in N_color:
    file = open("results/" + N + "_data_delta_E.txt", 'r')
    lines = file.readlines()
    lbl = "N = " + N
    X = []
    Y = []
    for i in range(7,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = c, label = lbl)



# N = "6"

# print("plotting Delta E ...")
# file = open("results/" + N + "_data_delta_E.txt", 'r')
# lines = file.readlines()
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(7,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'red', label = lbl)
# subfig1.set_xlabel(r'$J_1$ / $J_2$', fontsize = 18)
# subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 18)
# subfig1.set_title(r'Anregungsenergieren $\Delta E$', fontsize = 18)
# subfig1.axhline(0, color = "grey")

# N = "8"

# file = open("results/" + N + "_data_delta_E.txt", 'r')
# lines = file.readlines()
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(7,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)

# N = "10"

# file = open("results/" + N + "_data_delta_E.txt", 'r')
# lines = file.readlines()
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(7,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'green', label = lbl)


# N = "12"

# file = open("results/" + N + "_data_delta_E.txt", 'r')
# lines = file.readlines()
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(7,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'orange', label = lbl)

# N = "14"

# file = open("results/" + N + "_data_delta_E.txt", 'r')
# lines = file.readlines()
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(7,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'brown', label = lbl)

# N = "16"

# file = open("results/" + N + "_data_delta_E.txt", 'r')
# lines = file.readlines()
# lbl = "N = " + N
# X = []
# Y = []
# for i in range(7,len(lines)):
#     x, y = lines[i].split("\t")
#     #print(x + " " + y + "\r")
#     X += [float(x)]
#     Y += [float(y)]
# subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'purple', label = lbl)

subfig1.axhline(0, color = "grey")
subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("results/" + "delta_E.png")
#plt.show()

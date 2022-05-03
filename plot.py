import matplotlib.pyplot as plt

def plot_delta_E():
    file = open("results/data_delta_E.txt", 'r')
    lines = file.readlines()
    lbl = ""
    for i in range(4):
        lbl += lines[i]
    X = []
    Y = []
    for i in range(7,len(lines)):
        x, y = lines[i].split("\t")
        #print(x + " " + y + "\r")
        X += [float(x)]
        Y += [float(y)]
    fig1, subfig1 = plt.subplots(1,1,figsize=(16,9))
    subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)
    subfig1.set_xlabel(r'J_1 / J_2', fontsize = 14)
    subfig1.set_ylabel(r'$\Delta E = E_1 - E_0$  in $J_2$', fontsize = 14)
    subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)
    plt.savefig("results/data_delta_E.png")
    plt.show()

if __name__ == "__main__":
    print("plotting delta E")
    plot_delta_E()
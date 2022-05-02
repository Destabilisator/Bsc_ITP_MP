import matplotlib.pyplot as plt

file = open("data.txt", 'r')

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
#markers, caps, bars = subfig1.errorbar(X, Y, YErr, lw = 1, fmt = "o", ls = "solid", markersize = 0, color = 'red', ecolor = "blue", capsize = 2.0 , label = lbl)

subfig1.plot(X, Y, lw = 1, ls = "solid", markersize = 2, marker = "o", color = 'blue', label = lbl)

# transparent error bars
#[bar.set_alpha(0.25) for bar in bars]
#[cap.set_alpha(0.25) for cap in caps]

subfig1.set_xlabel('J1/J2', fontsize = 14)
subfig1.set_ylabel('E1 - E0 in J2', fontsize = 14)

subfig1.legend(loc = 'best' ,frameon = False, fontsize = 14)

plt.savefig("data.png")

plt.show()

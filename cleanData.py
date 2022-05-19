import os
import sys

N = sys.argv[1]

print("removing old data ...")
for filename in os.listdir("results/3DData/" + N + "/C/"):
    if filename == "dummyFile.txt":
        continue
    os.remove("results/3DData/" + N + "/C/" + filename)
    os.remove("results/3DData/" + N + "/X/" + filename)
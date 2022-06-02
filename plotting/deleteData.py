import os
import sys
import shutil

N = sys.argv[1]

print("removing old data (" + N + ") ...")
for filename in os.listdir("results/3DData/" + N + "/C/"):
    if filename == "dummyFile.txt":
        continue
    try:
        shutil.rmtree("results/3DData/" + N + "/C/" + filename, ignore_errors=False, onerror=False)
    except:
        os.remove("results/3DData/" + N + "/C/" + filename)

for filename in os.listdir("results/3DData/" + N + "/X/"):
    if filename == "dummyFile.txt":
        continue
    try:
        shutil.rmtree("results/3DData/" + N + "/X/" + filename, ignore_errors=False, onerror=False)
    except:
        os.remove("results/3DData/" + N + "/X/" + filename)
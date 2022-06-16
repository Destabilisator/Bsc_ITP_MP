import os
import sys
import shutil

N = sys.argv[1]

which =  sys.argv[2]

if which == "3D":
    print("removing old 3D data (" + N + ") ...")
    for filename in os.listdir("results/3DData/" + N + "/C/"):
        if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
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
            os.remove("results/" + N + "/X/" + filename)
elif which == "SG":
    print("removing old spin gap data (" + N + ") ...")
    for filename in os.listdir("results/" + N + "/data/spin_gap_data/"):
        if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
            continue
        try:
            shutil.rmtree("results/" + N + "/data/spin_gap_data/" + filename, ignore_errors=False, onerror=False)
        except:
            os.remove("results/" + N + "/data/spin_gap_data/" + filename)
elif which == "EE":
    print("removing old spin gap data (" + N + ") ...")
    for filename in os.listdir("results/" + N + "/data/excitation_energies_data/"):
        if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
            continue
        try:
            shutil.rmtree("results/" + N + "/data/excitation_energies_data/" + filename, ignore_errors=False, onerror=False)
        except:
            os.remove("results/" + N + "/data/excitation_energies_data/" + filename)
else:
    print("no valid data to delete")
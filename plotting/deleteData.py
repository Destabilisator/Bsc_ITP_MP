import os
import sys
import shutil

N = sys.argv[1]

which = sys.argv[2]

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
    for i in range(1,6):
        for filename in os.listdir("results/" + N + "/data/spin_gap_data/" + str(i) + "/"):
            if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
                continue
            try:
                shutil.rmtree("results/" + N + "/data/spin_gap_data/" + str(i) + "/" + filename, ignore_errors=False, onerror=False)
            except:
                os.remove("results/" + N + "/data/spin_gap_data/" + str(i) + "/" + filename)
elif which == "EE":
    print("removing old excitation energy data (" + N + ") ...")
    for i in range(1,6):
        for filename in os.listdir("results/" + N + "/data/excitation_energies_data/" + str(i) + "/"):
            if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
                continue
            try:
                shutil.rmtree("results/" + N + "/data/excitation_energies_data/" + str(i) + "/" + filename, ignore_errors=False, onerror=False)
            except:
                os.remove("results/" + N + "/data/excitation_energies_data/" + str(i) + "/" + filename)
elif which == "BRT":
    print("removing old benchmarking run time ...")
    for i in range(1,6):
        for filename in os.listdir("results/benchmarking/runtime/data/"):
            if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
                continue
            try:
                shutil.rmtree("results/benchmarking/runtime/data/" + filename, ignore_errors=False, onerror=False)
            except:
                os.remove("results/benchmarking/runtime/data/" + filename)
elif which == "BMU":
    print("removing old benchmarking memory usage ...")
    for i in range(1,6):
        for filename in os.listdir("results/benchmarking/memoryusage/data/"):
            if filename == "dummyFile.txt" or filename == "data_placeholder.txt":
                continue
            try:
                shutil.rmtree("results/benchmarking/memoryusage/data/" + filename, ignore_errors=False, onerror=False)
            except:
                os.remove("results/benchmarking/memoryusage/data/" + filename)
else:
    print("no valid data to delete")
    
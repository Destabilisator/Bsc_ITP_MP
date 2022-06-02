import sys
import glob
from PIL import Image
from typing import List
# import string

N = sys.argv[1]

print("creating gif (" + N + ", specific heat) ...")
rawName = N + "_specific_heat_"

def sort_filenames(fileList):
    length = len(fileList)
    for i in range(length-1):
        for j in range(0, length-i-1):
            name1 = fileList[j][len("./results/3DData/" + N + "_specific_heat_"): - len(".png")]
            name2 = fileList[j+1][len("./results/3DData/" + N + "_specific_heat_"): - len(".png")]
            if float(name1) > float(name2):
                fileList[j], fileList[j+1] = fileList[j+1], fileList[j]
    return fileList

frames = (Image.open(filename) for filename in sort_filenames(glob.glob("./results/3DData/" + N + "_specific_heat_*.png")))

frame = next(frames)

frame.save(fp = "results/3DData/" + rawName.rstrip("_") + ".gif", format = 'GIF', append_images = frames, save_all = True,
           duration = 200, loop = 0)

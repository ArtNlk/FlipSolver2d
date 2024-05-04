from PIL import Image
import random
import re
from collections import defaultdict
import sys

parseRegex = re.compile("(\d+):(-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*)")
    
def defaultValue():
    return None

inputFile = "data_eigen.txt"
width = 64
height = 64

def linIdxOfOffset(linidx, iOffset, jOffset):
    return linidx + iOffset * width + jOffset

def index2d(linIdx):
    output = []
    output[0] = linIdx // width
    output[1] = linIdx - output[0]*width
    return output

data = defaultdict(defaultValue)

picSize = width*height

# Create a new image with RGB mode
image = Image.new("RGB", (picSize, picSize))

# Create a pixel access object
pixels = image.load()

with open(inputFile) as f:
    for line in f:
        if "=" in line:
            continue
        vals = line.split(' ')
        row = int(vals[0])
        col = int(vals[1])
        val = float(vals[2])
        pixels[row,col] = (255,255,255)

#for x in range(picSize):
#    for y in range(picSize):
#        if pixels[x,y] != pixels[y,x]:
#            print("Asymmetric!")
# Save the image
image.save("viz_eigen.png")
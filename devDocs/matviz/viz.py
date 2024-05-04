from PIL import Image
import random
import re
from collections import defaultdict
import sys

parseRegex = re.compile("(\d+):(-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*)")
    
def defaultValue():
    return None

inputFile = "data.txt"
width = 0
height = 0

class Unit:
    index = 0
    diag = 0.0
    iPos = 0.0
    iNeg = 0.0
    jPos = 0.0
    jNeg = 0.0

def linIdxOfOffset(linidx, iOffset, jOffset):
    return linidx + iOffset * width + jOffset

def index2d(linIdx):
    output = []
    output[0] = linIdx // width
    output[1] = linIdx - output[0]*width
    return output

data = defaultdict(defaultValue)

with open(inputFile) as f:
    for line in f:
        if "=" in line:
            continue
        m = parseRegex.search(line)
        if not m:
            sizeList = line.split(',')
            if len(sizeList) == 2:
                height = int(sizeList[0])
                width = int(sizeList[1])
                print(f"Found size (WxH): {width} x {height}")
                continue
            else:
                print(f"Bad match: \"{line}\"")
                continue
        u = Unit()
        u.index = int(m.group(1))
        u.diag = float(m.group(2))
        u.iNeg = float(m.group(3))
        u.iPos = float(m.group(4))
        u.jNeg = float(m.group(5))
        u.jPos = float(m.group(6))
        data[u.index] = u


picSize = width*height

# Create a new image with RGB mode
image = Image.new("RGB", (picSize, picSize))

# Create a pixel access object
pixels = image.load()

for i in range(picSize):
    u = data[i]
    if not u:
        continue
    
    iPosIdx = linIdxOfOffset(i,1,0)
    iNegIdx = linIdxOfOffset(i,-1,0)
    jPosIdx = linIdxOfOffset(i,0,1)
    jNegIdx = linIdxOfOffset(i,0,-1)
    
    if i == 2052:
        pixels[i,0] = (255,0,0)
    
    pixels[i, i] = (255,255,255) if u.diag > 0.0 else (0,0,0)
    if 0 <= iPosIdx < picSize: pixels[i,iPosIdx] = (255,255,255) if abs(u.iPos) > 0.0 else (0,0,0)
    if 0 <= iNegIdx < picSize: pixels[i,iNegIdx] = (255,255,255) if abs(u.iNeg) > 0.0 else (0,0,0)
    if 0 <= jPosIdx < picSize: pixels[i,jPosIdx] = (255,255,255) if abs(u.jPos) > 0.0 else (0,0,0)
    if 0 <= jNegIdx < picSize: pixels[i,jNegIdx] = (255,255,255) if abs(u.jNeg) > 0.0 else (0,0,0)

#for x in range(picSize):
#    for y in range(picSize):
#        if pixels[x,y] != pixels[y,x]:
#            print("Asymmetric!")
# Save the image
image.save("viz.png")
from PIL import Image
import random
import re
from collections import defaultdict
import sys
from scipy.sparse import dok_matrix
import numpy as np

parseRegex = re.compile("(\d+):(-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*) (-?\d+\.*\d*)")
    
def defaultValue():
    return None

inputFile = "data.txt"
vecInFile = "vin.txt"
vecOutFile = "vout.txt"

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

vin = sparse_vector = np.zeros(width*height, dtype=np.float64)
vout = sparse_vector = np.zeros(width*height, dtype=np.float64)

with open(vecInFile) as f:
    i = 0
    for line in f:
        if "=" in line:
            continue
        vin[i] = float(line)
        i+=1

with open(vecOutFile) as f:
    i = 0
    for line in f:
        if "=" in line:
            continue
        vout[i] = float(line)
        i+=1



picSize = width*height

mat = dok_matrix((picSize, picSize), dtype=np.float64)

for i in range(picSize):
    u = data[i]
    if not u:
        continue
    
    iPosIdx = linIdxOfOffset(i,1,0)
    iNegIdx = linIdxOfOffset(i,-1,0)
    jPosIdx = linIdxOfOffset(i,0,1)
    jNegIdx = linIdxOfOffset(i,0,-1)
    
    mat[i, i] = u.diag
    if 0 <= iPosIdx < picSize: mat[i,iPosIdx] = u.iPos
    if 0 <= iNegIdx < picSize: mat[i,iNegIdx] = u.iNeg
    if 0 <= jPosIdx < picSize: mat[i,jPosIdx] = u.jPos
    if 0 <= jNegIdx < picSize: mat[i,jNegIdx] = u.jNeg

testVout = mat.dot(vin)

for i in range(picSize):
    if abs(testVout[i] - vout[i]) > 0.0001:
        print(f"Mismatch: {testVout[i]} vs {vout[i]}")
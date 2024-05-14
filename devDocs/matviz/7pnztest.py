from PIL import Image
import random
import re
from collections import defaultdict
import sys
from scipy.sparse import dok_matrix
import numpy as np

size = 64

iOffset = 1
jOffset = size
kOffset = size*size

picSize = size*size

mat = dok_matrix((picSize, picSize), dtype=np.float64)

for i in range(picSize):
    ipIdx = i + iOffset
    jpIdx = i + jOffset
    kpIdx = i + kOffset
    
    inIdx = i - iOffset
    jnIdx = i - jOffset
    knIdx = i - kOffset
    
    mat[i,i] = 1
    if 0 <= ipIdx < picSize: mat[i,ipIdx] = 1
    if 0 <= jpIdx < picSize: mat[i,jpIdx] = 1

out = mat * mat.transpose()

image = Image.new("RGB", (picSize, picSize))

# Create a pixel access object
pixels = image.load()

rows,cols = out.nonzero()
for row,col in zip(rows,cols):
    pixels[int(row),int(col)] = (255,255,255)

image.save("viz_ipp.png")
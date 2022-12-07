import numpy as np
import sys

inputfile = sys.argv[1]

with open(inputfile) as f:
   data =  f.readlines()

#print(data)
data = np.array(list(map(np.float64,data)))

print(np.min(data),np.max(data))

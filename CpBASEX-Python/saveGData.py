import numpy as np
import time
from pbasex import get_gData
import sys

gData = {}
gData['x'] = np.arange(2**9, dtype='double')
xkratio = int(sys.argv[1])
gData['l'] =  2*np.arange(2+1)
gData['params'] = 0.7*xkratio #sigma
gData['rBF'] = 'gauss'
gData['trapzStep'] = 0.1

if xkratio == 1:
	gData['k'] = gData['x']
	gData['k'][gData['k']==0] = gData['trapzStep']
else:
	gData['k'] = gData['x'][xkratio/2-1::xkratio] + 0.5

np.seterr("ignore")

t0 = time.time()
get_gData(gData)
print time.time()-t0

import numpy as np
import time
import sys
sys.path.append('../CpBASEX-Python')
from pbasex import get_gData

gData = {}
gData['x'] = np.arange(2**9, dtype='double') # 512x512 quadrant
xkratio = 1 # 64 Radial basis functions
lmax = 4 # Up to l=2
gData['l'] =  2*np.arange(lmax/2+1).astype(int)
gData['rBF'] = 'gauss' # Gaussian basis function
gData['params'] = 0.7*xkratio # Gaussian width
gData['trapzStep'] = 0.1 # Integration accuracy

if xkratio == 1:
	gData['k'] = gData['x']
	gData['k'][gData['k']==0] = gData['trapzStep']
else:
	gData['k'] = gData['x'][int(xkratio/2)-1::xkratio] + 0.5

np.seterr("ignore")

t0 = time.time()
get_gData(gData)
print(time.time()-t0)

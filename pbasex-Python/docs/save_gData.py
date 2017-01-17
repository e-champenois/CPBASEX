import numpy as np
import time
from pbasex import get_gData

# Settings
save_path = None # Automatic file naming if none
save_dir = None # Directory to save the data
nx = 512  # 512x512 quadrant
xkratio = 2 # Ratio of radial basis functions to pixel radii
lmax = 4 # Up to l=4
k_spacing = 'linear' # Even pixel or energy bins
gData = {}
gData['rBF'] = 'gauss' # Gaussian basis function

# Set up the parameters
gData['x'] = np.arange(nx, dtype='double')

if k_spacing=='linear':

	gData['k'] = np.arange(0, nx, xkratio) + 0.5 * (xkratio - 1)
	gData['params'] = 0.7 * xkratio # Gaussian width

elif k_spacing=='quadratic':

	gData['k'] = np.sqrt(np.linspace(0, (nx-1)**2, nx))
	gData['params'] = xkratio

gData['l'] =  2 * np.arange(lmax/2 + 1).astype(int)

if gData['rBF'] == 'custom':

	def rBF(r, k, params):

		return 1/((x-k)**2+(params/2)**2) # Lorentzian basis function

	def zIP(r, k, params):

		return np.sqrt((np.sqrt(max(0, 10-(params/2)**2)) + k)**2 - r**2)

	trapz_step = 0.05

	custom_rBF = (rBF, zIP, trapz_step)

else:

	custom_rBF = None

np.seterr("ignore")

t0 = time.time()
get_gData(gData, save_path=save_path, save_dir=save_dir, custom_rBF=custom_rBF)
print(time.time()-t0)
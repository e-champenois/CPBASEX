from __future__ import division
import numpy as np
import h5py
try:
	import dill
except:
	def cpu_count():
		return 1
else:
	from multiprocessing import Pool, cpu_count
	def packed_func(func_and_item):
		dumped_function, item = func_and_item
		target_function = dill.loads(dumped_function)
		res = target_function(item)
		return res
	def pack(target_function, items):
		dumped_function = dill.dumps(target_function)
		dumped_items = [(dumped_function, item) for item in items]
		return packed_func,dumped_items

def pbasex(images, gData, make_images=False, weights=None, regularization=0, alpha=1):
	"""Perform an Abel inversion using the pBASEX algorithm (see Garcia et al., Rev. Sci. Instrum. 75, 4989 (2004)).

	Parameters
	----------
	images : {2-3}D ndarray
		Pixelated folded images to Abel invert, stacked in the 3rd dimension.
	gData : str/h5py.File/dict
		Reference to the gData dictionary, holding the matrices needed for pBASEX inversion. See the loadG and get_gData functions for more information.
	make_images : bool
		Boolean describing whether or not to generate the fitted and inverted images. Defaults to False due to added computation.
	weights : 2D array
		Weights associated with the error at each image pixel for the least squares fitting. This adds some computation time to the algorithm.
	regularization : int/1D array
		Value(s) of the L-2 norm regularization parameter to use. Defaults to 0, indicating no regularization.
	alpha : int
		VMI spectrometer parameter correlating a radius to an energy, following E = alpha * r**2.

	Returns
	-------
	out : dict
		Dictionary holding the inversion outputs with key/value pairs:
		'E' : 1D array
			The energies sampled in the reconstructed spectrum.
		'IE' :{1-2}D ndarray
			The reconstructed kinetic energy spectra, stacked for different input images in the 2nd dimension.
		'betas' : {2-3}D ndarray
			The reconstructed beta values as a function of energy (1st dimension), l-order (2nd dimension), and input image (3rd dimension).
		'fit' {2-3D} ndarray
			The result of the input image fit. Only output if the make_images parameter is set to True.
		'inv' {2-3D} ndarray
			The inversion of the input image fit. Only output if the make_images parameter is set to True.

	See also
	--------
	foldQuadrant : Fold images by quadrant.
	loadG : Loads gData dictionary.

	Examples
	--------
	>>>images = imread(images_filename)
	>>>gData = loadG(gData_filename)
	>>>folded_images = foldQuadrant(images)
	>>>out = pbasex(folded_images, gData, make_images=True, regularization=1/gData['x'], alpha=4e-5)
	"""

	# Load the gData and find the dimension of the problem.
	gData = loadG(gData, make_images)
	nx, nk, nl = len(gData['x']), gData['nk'], gData['nl']
	try:
		nim = images.shape[2]
	except:
		nim = 1
		images = images.reshape(nx,nx,nim)

	# Invert the images through a least-squares fit of the Abel transformed basis functions.
	images = images.reshape(nx**2,nim)

	if weights is None:
		c = gData['V'].dot((gData['S']/(gData['S']**2+regularization))[:,None]*(gData['Up'].dot(images)))
	else:
		weights = weights.flatten()
		c = gData['V'].dot((gData['S']/(gData['S']**2+regularization))[:,None]*(np.linalg.solve((gData['Up']*weights[None,:]).dot(gData['Up'].T),gData['Up'].dot(weights[:,None]*images))))

	# Calculate kinetic energy spectra and angular distributions from the fit coefficients.
	E = alpha*gData['x']**2
	IEB = 1/(2*alpha)*np.diag(gData['x']).dot(gData['frk'].dot(c.reshape(nl,nk,nim).swapaxes(0,1).reshape(nk,nl*nim)))
	IE = IEB[:,:nim]
	betas = IEB[:,nim:].reshape(nx,nl-1,nim)/IE[:,None,:]

	if make_images:
		fit = unfoldQuadrant(gData['Up'].T.dot(np.diag((gData['S']**2+regularization)/gData['S']).dot(gData['V'].T.dot(c))).reshape(nx,nx,nim))
		inv = unfoldQuadrant(gData['Ginv'].dot(c).reshape(nx,nx,nim))

	# Populate the output dictionary
	out = {'E': E, 'IE': np.squeeze(IE), 'betas': np.squeeze(betas), 'c': np.squeeze(c)}
	if make_images:
		out['fit'], out['inv'] = np.squeeze(fit), np.squeeze(inv)

	return out

def loadG(gData, make_images=False):
	"""Loads the gData dictionary, holding the matrices needed for the pBASEX inversion.

	Parameters
	----------
	gData : str/h5py.File/dict
		loadG is called recursively from a string indicating the file path, to an open h5py.File, to a gData dictionary that is checked for format.
	make_images : bool
		Boolean describing whether or not to load the matrices needed to generate the fitted and inverted images.

	Returns
	-------
	gData : dict
		The gData dictionary used for pBASEX inversion, with key/value pairs:
		'x' : 1D array
			The pixel coordinates along one axis. A symmetric grid is assumed.
		'nk' : int
			The number of radial basis functions.
		'nl' : int
			The number of angular basis functions.
		'Up', 'S', 'V' : 2D ndarrays
			The matrices directly used for pBASEX inversion. See the get_gData function for more information.
		'frk' : 2D ndarray
			A matrix representing the sampled radial basis functions.

	See also
	--------
	pbasex : Perform an Abel inversion.
	get_gData : Calculate the gData dictionary.

	Examples
	--------
	>>> gData = loadG(gData_filename)
	>>> with hp5y.File(gData_filename, 'r') as h5file:
	>>>		gData = loadG(h5file)
	"""

	if isinstance(gData, str):
		with h5py.File(gData, 'r') as gData:
			return loadG(gData, make_images)
	elif isinstance(gData, h5py.File):
		gData_dict = {}
		for key in ['x','nk','nl','Up','S','V','frk']+make_images*['Ginv']:
			gData_dict[key] = gData[key].value
		return loadG(gData_dict, make_images)
	else:
		return gData

def get_gData(gData, h5_filename=None, custom_rBF=None, nProc=cpu_count()):
	"""Calculates and aggregates the data needed for the pBASEX algorithm, saving the results in an h5 file.

	Parameters
	----------
	gData : dict
		The base gData dictionary defining the basis functions, with key/value pairs:
			'x' : 1D array
				The pixel coordinates along one axis. A symmetric grid is assumed.
			'k' : 1D array
				A vector indexing the radial basis functions.
			'l' : 1D array
				A vector indexing the angular basis functions. This vector should start with the value 0 and contain only non-negative, even values.
			'params' : int/tuple/Python object
				The additional parameters passed to the radial basis functions.
			'rBF' : str
				Indicates which radial basis function to use. The keyword 'custom' will use the custom_rBF parameter as the radial basis function, while other strings will select the appropriate pre-defined function from the rBFs dictionary.
			'trapzStep' : float
				Indicates the step size used in the numerical integration for the Abel transformed basis functions.
	h5_filename : str
		Indicates a path to save the result to. Defaults to none, which saves an appropriately named file to the current path.
	custom_rBF : function/tuple
		Radial basis function to use in the pBASEX algorithm. This function must take three arguments, notably x (a position), k (indexing the radial basis functions), and params (other function paramters). This is an optional argument read used only when gData['rBF'] == 'custom'. This argument can also be a tuple of functions, where the first function will be used as the radial basis function and the second as the zero-integrand point function.
	nProc : int
		Optional integer indicating how many threads to parallelize the calculation to.

	See also
	--------
	findG : Numerically samples the Abel transformed basis functions.
	findGinv : Numerically samples the basis functions.

	Examples
	--------
	>>> gData = {'x':np.aranage(2**9, dtype='double'), 'k':np.arange(0.5,2**9,2), 'l':[0,2,4], 'params':1.4, 'rBF':'gauss', 'trapzStep':0.1}
	>>> get_gData(gData)
	"""

	# Set defaults.
	if h5_filename is None:
		h5_filename = 'G_r'+str(len(gData['x']))+'_k'+str(len(gData['k']))+'_l'+str(max(gData['l']))+'.h5'
	nProc = min(nProc, cpu_count())

	# Find problem dimension.
	gData['nk'] = len(gData['k'])
	gData['nl'] = len(gData['l'])

	# Load radial basis function and zero-integrand point function.
	if gData['rBF'] is not 'custom':
		rBF = rBFs[gData['rBF']]
	if callable(rBF) or len(rBF) == 1:
		zIP = 1.5*max(R)
		trapz_step = 0.05
	elif len(rBF) == 2:
		if callable(rBF[1]):
			rBF, zIP = rBF
			trapz_step = 0.05
		else:
			rBF, trapz_step = rBF
			zIP = 1.5*max(R)
	else:
		rBF, zIP, trapz_step = rBF
	rBF_2 = lambda x, k: rBF(x, k, gData['params'])
	zIP_2 = lambda x, k: zIP(x, k, gData['params'])

	# Sample the Abel transformed basis functions.
	G = findG(gData['x'], gData['k'], gData['l'], rBF_2, zIP_2, trapz_step, nProc)

	# Find the singular value decomposition of the G matrix for least-squares fitting.
	U, S, Vp = np.linalg.svd(G,0)
	gData['Up'], gData['S'], gData['V'] = U.T, S, Vp.T

	# Sample the radial basis functions.
	K, X = np.meshgrid(gData['k'], gData['x'])
	gData['frk'] = rBF_2(X, K)

	# Sample the basis functions.
	gData['Ginv'] = findGinv(gData['x'], gData['k'], gData['l'], rBF_2)

	# Save the calculation results.
	with h5py.File(h5_filename,'w') as f:
		for key in gData.keys():
			f.create_dataset(key,data=gData[key])

def findG(X, K, L, rBF, zIP, trapz_step, nProc):
	"""Numerically samples the Abel transformed basis functions.

	Parameters
	----------
	X : 1D array
		The pixel coordinates along one axis. A symmetric grid is assumed.
	K : 1D array
		A vector indexing the radial basis functions.
	L : 1D array
		A vector indexing the angular basis functions. This vector should start with the value 0 and contain only non-negative, even values.
	rBF : function
		A function describing the radial basis functions, with signature val = rBF(x, k).
	zIP : function
		A function describing at what radius to terminate integration, with signature val = zIP(r, k).
	trapz_step : float
		Indicates the step size used in the numerical integration for the Abel transformed basis functions.
	nProc : int
		Integer indicating how many threads to parallelize the calculation to.

	Returns
	-------
	G : 2D ndarray
		A matrix holding the values of the basis functions sampled on the pixel grid.

	See also
	--------
	get_gData : Calls findG while calculating and aggregating the data needed for the pBASEX algorithm.
	pack : Helps parallelize the computation.
	"""

	Y, X = np.meshgrid(X, X)
	Y, X = Y.flatten(), X.flatten()
	R = np.sqrt(X**2+Y**2)
	u = np.arange(0, max(zIP(0, K)), trapz_step)

	def findG_sub(yr):

		y, r = yr
		G_sub = np.zeros(len(K)*len(L))

		cos_term = np.nan_to_num(y/np.sqrt(u**2+r**2))
		leg_terms = np.zeros((len(L), len(u)))
		for i, l in enumerate(L):
			leg_terms[i] = leg(l, cos_term)

		for i, k in enumerate(K):
			u_sub = u[u<=zIP(r, k)]
			if u_sub.size:
				rad_term = rBF(np.sqrt(u_sub**2+r**2), k)
				G_sub[i::len(K)] = np.trapz(rad_term*leg_terms[:,:len(u_sub)], x=u_sub, dx=trapz_step)

		return G_sub

	if nProc > 1:
		p = Pool(nProc)
		G = p.map(*pack(findG_sub, zip(Y, R)))
		p.close()
	else:
		G = map(findG_sub, zip(Y, R))
	return np.array(list(G))/(2*np.pi)

def findGinv(X, K, L, rBF):
	#reoptimize??

	Y, X = np.meshgrid(X, X)
	Y, X = Y.flatten(), X.flatten()
	R = np.sqrt(X**2+Y**2)

	L, K = np.meshgrid(L, K)
	L, K = L.T.flatten(), K.T.flatten()

	CosTh = Y/R
	#CosTh[np.isnan(CosTh)] = 1

	R, K = np.meshgrid(R, K)

	return np.nan_to_num(rBF(R, K)*leg2D(L,CosTh)).T

def leg(l,x):

	if l==0:
		return np.ones(len(x))
	elif l==2:
		return (3*x**2-1)/2
	elif l==4:
		x2 = x**2
		return ((35*x2-30)*x2+3)/8
	else:
		raise ValueError('l-value not implemented.')

def leg2D(L,X):

	ls, idxs = np.unique(L, return_inverse=True)

	tmp = np.zeros((len(ls),len(X)))

	for i,l in enumerate(ls):
		tmp[i,:] = leg(l,X)

	out = np.zeros((len(L),len(X)))

	for i,idx in enumerate(idxs):
		out[i,:] = tmp[idx,:]

	return out

def foldQuadrant(M, x0=None, y0=None, quadrant_filter=[1,1,1,1]):

	big_int = 99999999
	qsigns = np.array([[-1,1],[-1,-1],[1,-1],[1,1]])
	try:
		sy,sx,sz = M.shape
	except ValueError:
		sy,sx = M.shape
		sz = 1
	if x0 is None:
		x0 = int(sx/2)
	if y0 is None:
		y0 = int(sy/2)
	quadrant_filter = np.array(quadrant_filter)

	lx,ly = big_int, big_int
	if quadrant_filter[1] or quadrant_filter[2]:
		lx = min(lx, x0+1)
	if quadrant_filter[0] or quadrant_filter[3]:
		lx = min(lx, sx-x0)
	if quadrant_filter[0] or quadrant_filter[1]:
		ly = min(ly, y0+1)
	if quadrant_filter[2] or quadrant_filter[3]:
		ly = min(ly, sy-y0)
	if sz>1:
		Mout = np.zeros((ly,lx,sz))
	else:
		Mout = np.zeros((ly,lx))

	for i in range(4):
		if quadrant_filter[i]:
			xf, yf = x0 + qsigns[i,1]*lx, y0 + qsigns[i,0]*ly
			if xf == -1:
				xf = None
			if yf == -1:
				yf = None
			Mout += M[y0:yf:qsigns[i,0], x0:xf:qsigns[i,1]]

	return Mout

def unfoldQuadrant(M):

	return np.vstack((np.hstack((np.rot90(M[1:,1:],2),np.flipud(M[1:,:]))),np.hstack((np.fliplr(M[:,1:]),M))))

def resizeFolded(M, r_max):

	try:
		sy,sx,sz = M.shape
		sz = (sz,)
	except ValueError:
		sy,sx = M.shape
		sz = ()

	if sx>r_max:
		x1, x2 = r_max, 0
	else:
		x1, x2 = sx, r_max-sx
	if sy>r_max:
		y1, y2 = r_max, 0
	else:
		y1, y2 = sy, r_max-sy

	return np.vstack((np.hstack((M[:y1,:x1],np.zeros((y1,x2)+sz))),np.zeros((y2,r_max)+sz)))

def gauss_rBF(x, k, params):

	return np.exp(-(x-k)**2/(2*params**2))/k**2

def gauss_zIP(r, k, params):

	return np.sqrt((np.sqrt(2*10)*params+k)**2-r**2)

rBFs = {'gauss': (gauss_rBF, gauss_zIP)}













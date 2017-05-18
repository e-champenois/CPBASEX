import numpy as np
from h5py import File
from os.path import join
from .rBFs import rBFs

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

def loadG(gData, make_images=False):

	if isinstance(gData, str):
		with File(gData, 'r') as gData:
			return loadG(gData, make_images)
	elif isinstance(gData, File):
		gData_dict = {}
		for key in ['x','nk','nl','Up','S','V','frk']+make_images*['Ginv']:
			gData_dict[key] = gData[key].value
		return loadG(gData_dict, make_images)
	else:
		return gData

def get_gData(gData, save_path=None, save_dir=None, custom_rBF=None, nProc=cpu_count()):

	# Set defaults.
	if save_path is None:
		save_path = '_'.join(['G','r'+str(len(gData['x'])),'k'+str(len(gData['k'])),'l'+str(max(gData['l']))])+'.h5'
	if save_dir is not None:
		save_path = join(save_dir, save_path)
	nProc = min(nProc, cpu_count())

	# Find problem dimension.
	gData['nk'] = len(gData['k'])
	gData['nl'] = len(gData['l'])

	# Load radial basis function and zero-integrand point function.
	if gData['rBF'] is not 'custom':
		rBF = rBFs(gData['rBF'], custom_rBF)
	if len(rBF) == 1:
		rBF = rBF[0]
	if callable(rBF):
		zIP = 1.5*max(gData['x'])
		trapz_step = 0.05
	elif len(rBF) == 2:
		if callable(rBF[1]):
			rBF, zIP = rBF
			trapz_step = 0.05
		else:
			rBF, trapz_step = rBF
			zIP = 1.5*max(gData['x'])
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
	with File(save_path,'w') as f:
		for key in gData.keys():
			f.create_dataset(key,data=gData[key])

def findG(X, K, L, rBF, zIP, trapz_step, nProc=cpu_count()):

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
	elif l==6:
		x2 = x**2
		return (((231*x2-315)*x2+105)*x2-5)/16
	elif l==8:
		x2 = x**2
		return ((((6435*x2-12012)*x2+6930)*x2-1260)*x2+35)/128
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
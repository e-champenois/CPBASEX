import numpy as np

def rBFs(type='gauss', custom_rBF=None):

	if type=='gauss' or type=='gaussian':
		return rBF_gauss, zIP_gauss

	elif type=='box':
		return rBF_box, zIP_box

	elif type=='custom':
		return custom_rBF

def rBF_gauss(r, k, params):

	return np.exp(-(r - k)**2 / (2*params**2))

def zIP_gauss(r, k, params):

	return np.sqrt((np.sqrt(2*10)*params + k)**2 - r**2)

def rBF_box(r, k, params):

	return (np.abs(r-k)<=params).astype(float)

def zIP_box(r, k, params):

	return np.sqrt((params + 1 + k)**2 - r**2)
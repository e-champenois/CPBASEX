import numpy as np
from .gData import loadG
from quadrant import unfoldQuadrant

def pbasex(images, gData, make_images=False, weights=None, regularization=0, alpha=1):
	
	gData  = loadG(gData, make_images)
	nx, nk, nl = len(gData['x']), gData['nk'], gData['nl']

	try:
		nim = images.shape[2]
	except:
		nim = 1
		images = images.reshape(nx, nx, nim)

	images = images.reshape(nx**2,nim)

	if weights is None:
		c = gData['V'].dot((gData['S']/(gData['S']**2+regularization))[:,None]*(gData['Up'].dot(images)))

	else:
		weights = weights.flatten()
		c = gData['V'].dot((gData['S']/(gData['S']**2+regularization))[:,None]*(np.linalg.solve((gData['Up']*weights[None,:]).dot(gData['Up'].T),gData['Up'].dot(weights[:,None]*images))))

	E = alpha*gData['x']**2
	IEB = 1/(2*alpha)*np.diag(gData['x']).dot(gData['frk'].dot(c.reshape(nl,nk,nim).swapaxes(0,1).reshape(nk,nl*nim)))
	IE = IEB[:,:nim]
	with np.errstate(divide='ignore'):
		betas = IEB[:,nim:].reshape(nx,nl-1,nim)/IE[:,None,:]

	if make_images:
		fit = unfoldQuadrant(gData['Up'].T.dot(np.diag((gData['S']**2+regularization)/gData['S']).dot(gData['V'].T.dot(c))).reshape(nx,nx,nim))
		inv = unfoldQuadrant(gData['Ginv'].dot(c).reshape(nx,nx,nim))

	out = {'E': E, 'IE': np.squeeze(IE), 'betas': np.squeeze(betas), 'c': np.squeeze(c)}
	if make_images:
		out['fit'], out['inv'] = np.squeeze(fit), np.squeeze(inv)

	return out


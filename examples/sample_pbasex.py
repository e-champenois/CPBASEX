import numpy as np
import matplotlib.pyplot as plt
from pbasex import pbasex, loadG
from quadrant import foldQuadrant, resizeFolded

# Load images
def load_image(path):
	return np.fromfile(path).reshape(1023,1023)
samples = ['low_counts', 'high_counts', 'fancy', 'fancier']
raw = np.dstack((load_image(sample+'.bin') for sample in samples))

# Fold images into quadrants
x0, y0 = 511, 511
quadrant_filter = [1,1,1,1]
folded = resizeFolded(foldQuadrant(raw, x0, y0, quadrant_filter), 512)

# Load inversion data
gData = loadG('G_r512_k512_l4.h5', 1)

# Apply the pBASEX algorithm
out = pbasex(folded, gData, make_images=True, alpha=4.1e-5)

# Plot some results
plt.figure(figsize=(12,9))
for i, sample in enumerate(samples):
	plt.subplot(4,4,4*i+1)
	plt.imshow(raw[:,:,i])
	plt.xticks([])
	plt.yticks([])
	plt.ylabel(sample)
	clim = plt.gci().get_clim()
	plt.clim(0,clim[1])
	if i==0:
		plt.title('Raw Image')
	plt.subplot(4,4,4*i+2)
	plt.plot(out['E'], out['IE'][:,i], 'k')
	plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
	if i==3:
		plt.xlabel('Energy (eV)')
	plt.gca().twinx()
	plt.plot(out['E'], out['betas'][:,:,i], '.', markersize=5)
	if i==0:
		plt.text(-3, 3.5, 'counts per eV', size='small')
		plt.text(12, 3.5, 'beta', size='small')
		plt.text(3.5, 3.25, 'I(E), ', color='black', size='large')
		plt.text(6, 3.25, 'B2', color='blue', size='large')
		plt.text(7.5, 3.25, ', ', color='black', size='large')
		plt.text(8, 3.25, 'B4', color='green', size='large')
	plt.ylim(-1,3)
	plt.subplot(4,4,4*i+3)
	plt.imshow(out['fit'][:,:,i]/4)
	plt.xticks([])
	plt.yticks([])
	plt.clim(0,clim[1])
	if i==0:
		plt.title('Fitted Image')
	plt.subplot(4,4,4*i+4)
	plt.imshow(out['inv'][:,:,i]/4)
	plt.xticks([])
	plt.yticks([])
	plt.clim(0,clim[1]/10)
	if i==0:
		plt.title('Inverted Image')
plt.show()
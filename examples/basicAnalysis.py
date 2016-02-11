import numpy as np
import sys
sys.path.append('../CpBASEX-Python')
from pbasex import pbasex, loadG, foldQuadrant, resizeFolded
import matplotlib as mpl
mpl.use('TKAgg')
import matplotlib.pyplot as plt

gData = loadG('G_r512_k64_l4.h5',1)

data = np.dstack((np.fromfile('high_counts.bin').reshape(1023,1023),
	   			  np.fromfile('low_counts.bin').reshape(1023,1023)))

x0,y0 = 512,512
fold = resizeFolded(foldQuadrant(data,x0,y0),512)

out = pbasex(fold,gData,1)

plt.figure()
plt.subplot(2,1,1)
plt.plot(out['E'],out['IE']/out['IE'].mean(0)[None,:])
plt.xlabel('Energy (eV)')
plt.ylabel('Percent of Total Counts')
plt.title('Photoelectron Spectrum')
plt.subplot(2,1,2)
plt.plot(out['E'],out['betas'][:,:,0])
plt.xlabel('Energy (eV)')
plt.ylabel('Anisotropy Parameters')
plt.title('Photoelectron Angular Distribution')
plt.figure()
plt.subplot(1,3,1)
plt.imshow(data[:,:,0])
plt.title('Data')
plt.xticks([])
plt.yticks([])
plt.subplot(1,3,2)
plt.imshow(out['fit'][:,:,0])
plt.title('Fitted Data')
plt.xticks([])
plt.yticks([])
plt.subplot(1,3,3)
plt.imshow(out['inv'][:,:,0])
plt.title('Inverted Fit')
plt.xticks([])
plt.yticks([])
plt.clim([0,1e4])
plt.show()
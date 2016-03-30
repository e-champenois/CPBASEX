import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../CpBASEX-Python')
from pbasex import pbasex, loadG, foldQuadrant, resizeFolded

gData = loadG('../gData/G_r512_k64_l4.h5',1)

simple_low_counts= np.fromfile('low_counts.bin').reshape(1023,1023)
simple_high_counts = np.fromfile('high_counts.bin').reshape(1023,1023)
fancy_high_counts = np.fromfile('fancy.bin').reshape(1023,1023)

data = np.dstack((simple_high_counts, simple_low_counts, fancy_high_counts))

x0,y0 = 511,511
fold = resizeFolded(foldQuadrant(data,x0,y0),512)

out = pbasex(fold, gData, make_images=True, alpha=4e-5)

plt.figure()
plt.subplot(2,1,1)
plt.plot(out['E'],(out['IE']/out['IE'].mean(0)[None,:])[:,:2])
plt.xlabel('Energy (eV)')
plt.ylabel('Percent of Total Counts')
plt.title('Photoelectron Spectrum')
plt.subplot(2,1,2)
plt.plot(out['E'],out['betas'][:,:,0])
plt.xlabel('Energy (eV)')
plt.ylabel('Anisotropy Parameters')
plt.title('Photoelectron Angular Distribution')
plt.figure()
plt.subplot(2,1,1)
plt.plot(out['E'],(out['IE']/out['IE'].mean(0)[None,:])[:,2])
plt.xlabel('Energy (eV)')
plt.ylabel('Percent of Total Counts')
plt.title('Photoelectron Spectrum')
plt.subplot(2,1,2)
filt = (out['IE']/out['IE'].mean(0)[None,:])[:,2]>(out['IE']/out['IE'].mean(0)[None,:])[:,2].max()*0.05
plt.plot(out['E'][filt],out['betas'][filt,:,2])
plt.xlabel('Energy (eV)')
plt.ylabel('Anisotropy Parameters')
plt.title('Photoelectron Angular Distribution')
plt.figure()
plt.subplot(1,3,1)
plt.imshow(data[:,:,2])
plt.title('Data')
plt.xticks([])
plt.yticks([])
plt.subplot(1,3,2)
plt.imshow(out['fit'][:,:,2])
plt.title('Fitted Data')
plt.xticks([])
plt.yticks([])
plt.subplot(1,3,3)
plt.imshow(out['inv'][:,:,2])
plt.title('Inverted Fit')
plt.xticks([])
plt.yticks([])
plt.clim([0,5e2])
plt.show()
import numpy as np
import sys
sys.path.append('../CpBASEX-Python')
from pbasex import pbasex, loadG, foldQuadrant, resizeFolded
import matplotlib as mpl
mpl.use('TKAgg')
import matplotlib.pyplot as plt

gData = loadG('G_r512_k128_l4.h5',1)

data = np.dstack((np.fromfile('high_counts.bin').reshape(1023,1023),
	   			  np.fromfile('low_counts.bin').reshape(1023,1023)))

x0,y0 = 512,512
fold = resizeFolded(foldQuadrant(data,x0,y0),512)

out = pbasex(fold,gData,1)

plt.ion()
plt.figure()
plt.subplot(2,1,1)
plt.plot(out['E'],out['IE']/out['IE'].mean(0)[None,:])
plt.subplot(2,1,2)
plt.plot(out['E'],out['betas'][:,:,0])
plt.figure()
plt.subplot(1,3,1)
plt.imshow(data[:,:,0])
plt.subplot(1,3,2)
plt.imshow(out['fit'][:,:,0])
plt.subplot(1,3,3)
plt.imshow(out['inv'][:,:,0])
plt.clim([0,1e3])
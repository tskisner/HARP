
import sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


Afile = 'A.out'
A = sio.mmread(Afile)
Adense = A.todense()
Apng = 'A.png'
plt.matshow(Adense)
plt.savefig ( Apng )
plt.clf()

invCfile = 'invC.out'
invC = sio.mmread(invCfile)
invCdense = invC.todense()
invCpng = 'invC.png'
plt.matshow(invCdense)
plt.savefig ( invCpng )
plt.clf()

for pix in ['signal', 'noise', 'measured']:
    datafile = pix + '.out'
    data = sio.mmread(datafile)
    png = pix + '.png'
    plt.matshow ( data )
    plt.savefig ( png )
    plt.clf()

truthfile = 'truth.out'
truth = sio.mmread(truthfile)
tpng = 'truth.png'

indx = np.arange ( 0, truth.shape[0], 1 )
plt.plot( indx, truth, linewidth=1 )
plt.savefig ( tpng )
plt.clf()




import sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


#Afile = 'A.out'
#A = sio.mmread(Afile)
#Adense = A.todense()
#Apng = 'A.png'
#plt.matshow(Adense)
#plt.savefig ( Apng )
#plt.clf()

#invCfile = 'invC.out'
#invC = sio.mmread(invCfile)
#invCdense = invC.todense()
#invCpng = 'invC.png'
#plt.matshow(invCdense)
#plt.savefig ( invCpng )
#plt.clf()

min = {}
max = {}
min['signal'] = 0.0
max['signal'] = 50.0
min['noise'] = -10.0
max['noise'] = 10.0
min['data'] = min['signal']
max['data'] = max['signal']

for pix in ['signal', 'noise', 'data']:
    datafile = pix + '.out'
    data = sio.mmread(datafile)
    png = pix + '.png'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow( data, vmin=min[pix], vmax=max[pix] )
    ax.set_aspect('auto')
    #fig.colorbar( cax )
    fig.savefig( png )

truthfile = 'truth.out'
truth = sio.mmread(truthfile)
tpng = 'truth.png'

indx = np.arange ( 0, truth.shape[0], 1 )

fig = plt.figure( )
ax = fig.add_subplot(111)
ax.plot( indx, truth, linewidth=1 )
ax.set_aspect('auto')
fig.savefig ( tpng )



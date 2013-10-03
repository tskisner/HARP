import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

file = sys.argv[1]
bins = int(sys.argv[2])
npix = int(sys.argv[3])

psf = np.zeros( [ bins, npix ] )

prow, pcol, pval = np.loadtxt( file, usecols=(0,1,2), unpack=True )

for vals in zip ( prow, pcol, pval ):
    psf[vals[0],vals[1]] += vals[2]

plt.figure(1)

plt.imshow( psf, interpolation='nearest', cmap=cm.Oranges, aspect='auto', vmin=0, vmax=1.0e-6 )

mappng = file + '.png'
plt.savefig ( mappng )
plt.clf()

plt.close()

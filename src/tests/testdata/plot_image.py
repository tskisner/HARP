import sys
import pyfits
import numpy as np
import pylab as py
import matplotlib.cm as cm


file = sys.argv[1]

imgdat = pyfits.getdata(file, 0)
imgnse = pyfits.getdata(file, 1)

fig = py.figure( figsize=(100, 12), dpi=600 )

frame1 = fig.add_subplot(211)
py.title( "Pixel Data (PSF projected input spectra + noise)" )
im1 = frame1.imshow(imgdat, aspect='equal', cmap=cm.get_cmap('jet', 100), vmin=-10, vmax=50)
#fig.colorbar(im1)

frame2 = fig.add_subplot(212)
py.title( "Inverse Pixel Noise Variance" )
im2 = frame2.imshow(imgnse, aspect='equal', cmap=cm.get_cmap('jet', 100), vmin=-10, vmax=50)
#fig.colorbar(im2)

py.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

# frame2 = fig.add_subplot(212)
# im2 = frame2.imshow(hits, aspect='equal', cmap=cm.get_cmap('gist_heat', 100), vmin=0, vmax=500)
# fig.colorbar(im2)
# py.xlim(100,250)
# py.ylim(150,250)
# im2.axes.get_xaxis().set_visible(False)
# im2.axes.get_yaxis().set_visible(False)

out = file + '.png'
py.savefig(out)


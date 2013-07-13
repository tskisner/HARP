import sys
import pyfits
import numpy as np
import pylab as py
import matplotlib.cm as cm


sigfile = 'sandbox_signal.fits.out'
nsefile = 'sandbox_noise.fits.out'
msrfile = 'sandbox_measured.fits.out'

imgsig = pyfits.getdata(sigfile, 0)
imgnse = pyfits.getdata(nsefile, 0)
imgmsr = pyfits.getdata(msrfile, 0)

fig = py.figure( figsize=(50, 6), dpi=600 )

frame1 = fig.add_subplot(311)
py.title( "Pixel Signal (PSF projected input spectra)" )
im1 = frame1.imshow(imgsig, aspect='equal', cmap=cm.get_cmap('jet', 100), vmin=-10, vmax=50)
#fig.colorbar(im1)

frame2 = fig.add_subplot(312)
py.title( "Pixel Noise" )
im2 = frame2.imshow(imgnse, aspect='equal', cmap=cm.get_cmap('jet', 100), vmin=-10, vmax=50)
#fig.colorbar(im2)

frame3 = fig.add_subplot(313)
py.title( "Pixel Data" )
im3 = frame3.imshow(imgmsr, aspect='equal', cmap=cm.get_cmap('jet', 100), vmin=-10, vmax=50)
#fig.colorbar(im3)

py.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

# frame2 = fig.add_subplot(212)
# im2 = frame2.imshow(hits, aspect='equal', cmap=cm.get_cmap('gist_heat', 100), vmin=0, vmax=500)
# fig.colorbar(im2)
# py.xlim(100,250)
# py.ylim(150,250)
# im2.axes.get_xaxis().set_visible(False)
# im2.axes.get_yaxis().set_visible(False)

py.savefig('img_projected.png')


import sys
import pyfits
import numpy as np
import pylab as py
import matplotlib.cm as cm


file = 'sandbox_results.fits.out'

imgtruth = pyfits.getdata(file, 0)
imgrt = pyfits.getdata(file, 1)
imgrf = pyfits.getdata(file, 2)
imgerr = pyfits.getdata(file, 3)

fig = py.figure( figsize=(10, 24), dpi=300 )

frame1 = fig.add_subplot(131)
py.title( "Input Spectra" )
im1 = frame1.imshow(imgtruth, aspect='equal', cmap=cm.get_cmap('jet', 100) )
fig.colorbar(im1)

frame2 = fig.add_subplot(132)
py.title( "Resolution Convolved Input" )
im2 = frame2.imshow(imgrt, aspect='equal', cmap=cm.get_cmap('jet', 100) )
fig.colorbar(im2)

frame3 = fig.add_subplot(133)
py.title( "Solution" )
im3 = frame3.imshow(imgrf, aspect='equal', cmap=cm.get_cmap('jet', 100) )
fig.colorbar(im3)

# frame2 = fig.add_subplot(212)
# im2 = frame2.imshow(hits, aspect='equal', cmap=cm.get_cmap('gist_heat', 100), vmin=0, vmax=500)
# fig.colorbar(im2)
# py.xlim(100,250)
# py.ylim(150,250)
# im2.axes.get_xaxis().set_visible(False)
# im2.axes.get_yaxis().set_visible(False)

py.savefig('result_spectra.png')


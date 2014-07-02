
import re
#import numpy as np

class HarpTest ( object ):
	'''Class for testing embedded python.'''

	def __init__ ( self, path, minval, maxval ):

		self.path = path
		self.minval = minval
		self.maxval = maxval
		self.diffval = maxval - minval

	def get_path ( self ):
		return self.path

	def get_diff ( self ):
		return self.diffval



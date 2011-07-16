"""
Matrix utilities.
"""
import sys

def print_matrix(mat, ignoreList=False, ignoreChar='.'):
	"""
	Print a matrix.
	At present this is adapted to printing connection matrices and
	replaces Infinity with its unicode representation. It might be good
	to generalize this later.

	There is also a "value ignore"/replace functionality, but it is not
	well developed at present. Matching or regex would be better, and a
	dictionary of 'find'=>'replacement' would also be nice.
	"""

	# TODO: Won't print matrices with a dimension beyond 99.
	# (Not that the terminal is ideal for that anyway.)

	# If there are any values we shouldn't print. 
	# XXX: This functionality is not very extensive. 
	doIgnore = True 
	if type(ignoreList) == type(False): # XXX: I should know better than this..
		doIgnore = False
	elif type(ignoreList) != list:
		ignoreList = [ignoreList]
	
	length = len(mat)
	width = len(mat[0])
	inf = float('Infinity')
	uinf = u'\u221e'

	# Header column numbers
	ln = " "*3
	if width < 10:
		ln = " "*2
	for i in range(width):
		if i < 10 or i %2 == 0:
			ln += "%d " % i
		else:
			ln += " "
	print ln

	# Matrix data
	for i in range(length):
		# Row headers. 
		if width < 10 or i >= 10:
			ln = "%d " % i
		else:
			ln = "%d  " % i

		# Row values
		for j in range(width):
			val = mat[i][j]
			if doIgnore and val in ignoreList:
				val = ignoreChar

			if val in [inf, -inf]:
				val = uinf

			if type(val) == bool:
				val = int(val)

			if type(val) not in [str, unicode]:
				val = str(val)

			ln += val + " "

		print ln


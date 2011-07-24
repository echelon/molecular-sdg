"""
Python support for enums.
"""

def enum(*sequential, **named):
	"""
	Python support for enums.

	Usage:
		>>> Numbers = enum('ZERO', 'ONE', 'TWO')
		>>> Numbers.ZERO
		0
		>>> Numbers.ONE
		1

	From: http://stackoverflow.com/questions/36932/whats-the-best-way-t
	o-implement-an-enum-in-python/1695250#1695250
	"""
	enums = dict(zip(sequential, range(len(sequential))), **named)
	return type('Enum', (), enums)

"""
Chain-related data structures and functions.
"""

class Chain(tuple):
	"""
	Chain Data Structure
	Contains chain atom path, metadata, etc.
	"""
	def __new__(cls, path):
		"""
		Build tuple subclass instance. The underlying tuple holds the
		chain path itself, so we do some cleanup to ensure each atom
		only occurs once. (Technically start or end could be the head.)
		"""
		atomsFound = {}
		for atom in path:
			if atom in atomsFound:
				raise Exception, "Chains cannot have duplicated atoms."
			atomsFound[atom] = True

		return tuple.__new__(cls, path)

	def __init__(self, path):
		"""
		Chain Constructor
		Must specify the atom path for the chain. 
		"""
		self.caps = [-1, -1] # There are two end caps. 
		self.zigzag = [] # Zigzag pattern along the chain from {L,R}.
		self.invertOk = False # Whether the chain may be inverted.  

	def __repr__(self):
		return "chain[%d, %s, %d]" % (self.caps[0], str(self[:]), self.caps[1])

	def __str__(self):
		return repr(self)

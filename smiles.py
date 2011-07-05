# FIXME: I don't have time to fix circular imports now
#from matrix import *

class Smiles(object):
	"""
	SMILES (Simplified Molecular Input Line Entry Specification)

	This class parses SMILES strings into a collection of tokens that
	are easily converted into molecular graphs.
	"""

	def __init__(self, smilesStr):
		"""CTOR."""

		self.string = smilesStr
		self.tokens = self.tokenizeString(smilesStr)

		# TODO: Canonoicalization
		self._cachedCanonical = None # Should be 'Smiles' instance

	def numAtoms(self):
		"""Calculates the number of non-H atoms"""
		cnt = 0
		for x in self.tokens:
			x = str(x)
			if x.isalpha() and x.upper() != 'h':
				cnt += 1
		return cnt

	def __str__(self):
		"""String Representation."""
		return "<Smiles: %s>" % self.string

	@staticmethod
	def tokenizeString(smileStr):
		"""
		Convert SMILES string into a list of tokens. 
		Important because 'Cl', 'Br', etc. constitute one token.
		"""
		# TODO: This won't work for inorganic set (bracketed) 
		# TODO: Stereochemistry
		# TODO: This won't work for ring closures beyond the 9th.
		READ_AHEAD = {
			'C': 'Cl',
			'B': 'Br'
		}

		queue = []
		pos = 0
		while pos < len(smileStr):
			ch = smileStr[pos]

			# Check for Br, Cl, etc.
			if ch in READ_AHEAD and pos < len(smileStr)-1:
				rd = ch + smileStr[pos+1]
				if rd == READ_AHEAD[ch]:
					ch = rd
					pos += 1

			# TODO: Connectivity beyond '9'

			queue.append(ch)
			pos += 1

		return queue

	@staticmethod
	def canonicalize(smiles):
		"""
		Adapted from Morgan's original algorithm as presented in
		[TODO: Source.]

		Input: a Smiles() object.
		Output: a Smiles() object.
		"""

		def computeInvariant():
			"""
			"""
			pass


		# We need graph invariants (topological indices not dependant 
		# on the labeling scheme) and labels for each atom.
		invariants = [1 for x in range(len(smiles.numAtoms()))]
		labels = [0 for x in range(len(smiles.numAtoms()))]


		invar = None# TODO: Compute invariant
		l = []


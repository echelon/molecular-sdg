from matrix import *

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
		# XXX: Perhaps this isn't something I should even manage here...
		self._cachedCanonical = None # Should be 'Smiles' instance

	def numAtoms(self):
		"""Calculates the number of non-H atoms"""
		# XXX: This hydrogen-excluded policy WILL cause problems in the 
		# future if I don't develop a comprehensive protocol for
		# handling.
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

	def toMatrix(self):
		"""Convert a SMILES string into a adjacency matrix representation."""

		# First, we must convert the input string into a proper queue
		# Organic subset: B, C, N, O, P, S, F, Cl, Br, I 
		# Everything else must be specified in brackets. 
		ATOMS = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
		ATOMS_UPPER = map(lambda x: x.upper(), ATOMS)

		# Additional Symbols
		BRANCHING = '()' # Branch Start & End
		BONDS = '=#'	 # Double and Triple Bonds
		CONNECTIVE = '%' # Connectivity beyond '9' -- TODO NOT YET HANDLED.

		queue = self.tokens
		mat = MolMatrix(self.numAtoms())

		# Symbol type tests
		def is_atom(sym): return sym.upper() in ATOMS_UPPER
		def is_bond(sym): return sym in BONDS
		def is_branch_start(sym): return sym == '('
		def is_branch_end(sym): return sym == ')'
		def is_closure(sym): return sym.isdigit()
		def is_bond_order(sym): return sym in '=#'

		# Make note of the connection beween two atoms. 
		def connect(a1, a2, bondOrder=1):
			mat.connectMat[a1][a2] = bondOrder
			mat.connectMat[a2][a1] = bondOrder

		# Keep track of state during parsing
		atomCnt = 0 # Total num atoms
		atomPrev = 0 # Previous atom 
		bondOrder = 1 # Bond order: 1, 2, 3
		branchStack = [] # Branching, eg. C(C)(C)C
		ringClosures = {} # Keep track of cycles, eg. c1ccccc1

		for sym in queue:

			if is_atom(sym):
				atomId = atomCnt
				atomCnt += 1

				if atomId != 0:
					connect(atomPrev, atomId, bondOrder)
					bondOrder = 1

				atomPrev = atomId

			elif is_branch_start(sym):
				branchStack.append(atomPrev)

			elif is_branch_end(sym): 
				atomPrev = branchStack.pop()

			elif is_closure(sym): # TODO: Won't handle > 9
				num = int(sym)
				if num not in ringClosures:
					ringClosures[num] = atomPrev
				else:
					connect(atomPrev, ringClosures[num], bondOrder)
					bondOrder = 1

			elif is_bond_order(sym):
				if sym == '=':
					bondOrder = 2
				else:
					bondOrder = 3

		# End queue processing, return matrix.
		return mat

def smiles_to_matrix(smiles):
	"""Function to return a MolMatrix from a smiles string without an
	intermediate."""
	# XXX: Do I really need to keep this?

	if type(smiles) == str:
		smiles = Smiles(smiles)

	return smiles.toMatrix()


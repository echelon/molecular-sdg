#!/usr/bin/env python

from smiles import Smiles

# XXX: Tis will be messy until I finish porting the algorithm
def smiles_to_matrix(string):
	"""Convert a SMILES string into a adjacency matrix representation."""

	# First, we must convert the input string into a proper queue
	# Organic subset: B, C, N, O, P, S, F, Cl, Br, I 
	# Everything else must be specified in brackets. 
	ATOMS = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
	ATOMS_UPPER = map(lambda x: x.upper(), ATOMS)

	BRANCHING = '()' # Branch Start & End
	BONDS = '=#'	 # Double and Triple Bonds
	CONNECTIVE = '%' # Connectivity beyond '9' -- TODO NOT YET HANDLED.

	smiles = Smiles(string)
	queue = smiles.tokens
	molmatrix = MolMatrix(smiles)
	matrix = molmatrix.matrix

	# Symbol type tests
	def is_atom(sym): return sym.upper() in ATOMS_UPPER
	def is_bond(sym): return sym in BONDS
	def is_branch_start(sym): return sym == '('
	def is_branch_end(sym): return sym == ')'
	def is_closure(sym): return sym.isdigit()
	def is_bond_order(sym): return sym in '=#'

	# Make note of the connection beween two atoms. 
	def connect(a1, a2, bondOrder=1):
		matrix[a1][a2] = bondOrder
		matrix[a2][a1] = bondOrder

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

	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE

	molmatrix.print_matrix()

class MolMatrix(object):
	"""Adjacency Matrix for molecules."""

	def __init__(self, molInput):
		"""CTOR"""

		if type(molInput) == str:
			molInput = Smiles(molInput)

		self.smiles = molInput
		
		sz = molInput.numAtoms()

		self.matrix= [[False for x in range(sz)] for xx in range(sz)]


	def print_matrix(self):
		"""Print the matrix. Debug."""
		print "AdjMat for %s" % self.smiles
		# Won't print > 100 atoms nicely
		ln = " "*3 if self.smiles.numAtoms() < 10 else " "*4

		# Header numbers
		for i in range(len(self.matrix)):
			if i < 10 or i %2 == 0:
				ln += "%d " % i
			else:
				ln += " "
		print ln

		# Graph data
		for i in range(len(self.matrix)):
			if self.smiles.numAtoms() < 10 or i >= 10:
				ln = "%d  " % i
			else:
				ln = "%d   " % i 
			for j in range(len(self.matrix)):
				ln += str(self.matrix[i][j]) + " " if self.matrix[i][j] else ". "
			print ln

# Temporary, for testing!!!
if __name__ == '__main__':
	import sys

	if len(sys.argv) < 2:
		print "Need to supply SMILES text as argument."
		sys.exit()
		
	mol = MolMatrix(sys.argv[1])
	smiles_to_matrix(sys.argv[1])

	#mol.print_matrix()	

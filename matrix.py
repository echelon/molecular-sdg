#!/usr/bin/env python

from smiles import Smiles

# XXX: This will be messy until I finish porting the algorithm
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

	return molmatrix

class MolMatrix(object):
	"""Adjacency Matrix for molecules."""

	def __init__(self, molInput):
		"""CTOR"""

		if type(molInput) == str:
			molInput = Smiles(molInput)

		self.smiles = molInput
		
		sz = molInput.numAtoms()

		self.matrix = [[False for x in range(sz)] for xx in range(sz)]

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

	def numAtoms(self):
		"""Reports the number of atoms in the molecule."""
		return len(self.matrix)

	def canonicalize(self):
		"""
		Adapted from Morgan's original algorithm as presented in
		[TODO: Source.]

		1. Label each node with its degree.
		2. Count number of different values.
		3. Recalculate labels by summing label values at neighbor nodes
		4. Count number of different values.
		5. Repeat from 3 until there is no change in number of different vals
			Awesome!!

		"""

		def num_distinct(li):
			"""Get the number of distinct vals in a list."""
			m = {}
			for x in li:
				m[x] = True
			return len(m)

		def get_neighbors(mat, atom):
			"""Get the neighbors of atom i."""
			n = []
			for i in range(len(mat)):
				if mat[atom][i]:
					n.append(i)
			return n

		def max_pos(li):
			"""Find the maximum position in a list"""
			hi = li[0]
			hiPos = 0
			for i in range(len(li)):
				if li[i] > hi:
					hi = li[i]
					hiPos = i
				elif li[i] == hi:
					pass # TODO: Bond orders to differentiate?
			return hiPos

		matrix = self.matrix

		size = len(matrix)
		newMatrix = [[False for x in range(size)] for xx in range(size)]

		# Atom labels -- used to build the new indices
		labels = [0 for x in range(size)]

		# Initially label each atom with its degree
		for i in range(size):
			for j in range(size):
				if matrix[i][j]:
					labels[i] += 1

		# Calculate atom labels until the number of distinct labels no
		# longer changes
		diff = num_distinct(labels)
		diff2 = diff
		repeatLoop = 5
		iterCount = 0
		while repeatLoop > 0 and iterCount < 1000:
			newLabels = [0 for x in range(size)]
			# Recompute label as the sum of the neighbor atom's labels.
			for i in range(size):
				for j in range(size):
					if not matrix[i][j]: 
						continue
					newLabels[i] += labels[j]

			labels = newLabels
			diff2 = num_distinct(labels)

			# XXX: Hack -- Continue at least 5 times after the values
			# appear to have converged.
			if diff == diff2:
				repeatLoop -= 1
				#break
			diff = diff2
			iterCount += 1


		# Score / Order the atoms by their labels.
		canonicalOrder = [] #[0 for x in range(size)]

		# Node One will be the maximum labelled atom.
		queue = []
		queue.append(max_pos(labels))

		# Canonicalization procedure. 
		while len(queue) > 0:
			front = queue.pop(0)
			canonicalOrder.append(front)
			nbrs = get_neighbors(matrix, front)
			while len(nbrs) > 0:
				# TODO: Need to use bond orders!
				hi = nbrs.pop(max_pos(nbrs))
				if hi not in canonicalOrder \
					and hi not in queue:
						queue.append(hi)

		# Build new matrix
		for i in range(size):
			newAtomPos = i
			oldAtomPos = canonicalOrder[i]
			neighbors = get_neighbors(matrix, oldAtomPos)

			for n in neighbors:
				oldNeighborPos = n
				newNeighborPos = canonicalOrder.index(n)

				newMatrix[newAtomPos][newNeighborPos] = \
						matrix[oldAtomPos][oldNeighborPos]

				newMatrix[newNeighborPos][newAtomPos] = \
						matrix[oldNeighborPos][oldAtomPos]


		# TODO TODO -- Don't replace current object!
		self.matrix = newMatrix

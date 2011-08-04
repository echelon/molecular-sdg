from molecule import Molecule
from atom import Atom

class MolMatrix(object):
	"""Adjacency Matrix for molecules."""

	def __init__(self, size):
		"""CTOR"""

		# Number of atoms.
		self.size = size

		# Connectivity matrix 
		self.connectMat = [[False for x in range(size)]
								for xx in range(size)]

		# Bond orders between atom pairs
		# 1/Single '-' and implicit, 2/Double '=', 3/Triple '#', 
		# TODO Aromatic as 1.5
		self.bondOrderMat = [[False for x in range(size)] 
									for xx in range(size)]

		# Atom Labels. 
		self.atomTypes = [False for x in range(size)] # C, O, N, Cl, etc.
		self.atomCharges = [0 for x in range(size)]
		self.atomIsotopes = [0 for x in range(size)]

	def getConnectMat(self):
		"""
		Build a connection matrix we can use in Floyd's algorithm.
		'Not connected' is represented by infinity.
		"""
		# TODO: Move to module where this is actually used...

		inf = float('Infinity')
		newMat = [[inf for x in range(self.size)] for xx in range(self.size)]

		for i in range(self.size):
			for j in range(self.size):
				val = self.connectMat[i][j]
				val = inf if not val else 1
				newMat[i][j] = val

		return newMat

	def print_matrix(self):
		"""
		Print the molecular data for debugging.
		"""
		print "Types: %s" % str(self.atomTypes)
		print "Charges: %s" % str(self.atomCharges)
		print "Isotopes: %s" % str(self.atomIsotopes)

		self.getHybridization(0) # XXX: Caches hybridizations. 
		print "Hybridizations: %s" % str(self._atomHybridizations)

		# XXX: Caches degrees
		for n in range(self.size):
			self.getDegree(n)

		print "Degrees: %s" % str(self._atomDegrees)
		print ""

		#print "AdjMat for %s" % self.smiles
		# XXX: Won't print >= 100 atoms nicely. Not that I would want
		# to print out such systems in the terminal...
		ln = " "*5 if self.size < 10 else " "*6

		# Header atoms
		for i in range(self.size):
			if len(self.atomTypes[i]) > 1:
				ln += self.atomTypes[i]
			else:
				ln += "%s " % self.atomTypes[i]

		print ln
		ln = " "*5 if self.size < 10 else " "*6

		# Header numbers
		for i in range(len(self.connectMat)):
			if i < 10 or i %2 == 0:
				ln += "%d " % i
			else:
				ln += " "
		print ln

		# Graph data
		for i in range(len(self.connectMat)):
			atom = self.atomTypes[i] # XXX: Can't print 3 char atoms...
			if len(atom) < 2:
				atom += " "
			if self.size < 10 or i >= 10:
				ln = "%s %d  " % (atom, i)
			else:
				ln = "%s %d   " % (atom, i) 
			for j in range(len(self.connectMat)):
				ln += str(int(self.bondOrderMat[i][j])) + " " \
						if self.bondOrderMat[i][j] else ". "
			print ln

	def canonicalize(self):
		"""
		Adapted from Morgan's original algorithm as presented in
		[Handbook of Cheminformatics Algorithms]

		1. Label each node with its degree.
		2. Count number of different values.
		3. Recalculate labels by summing label values at neighbor nodes
		4. Count number of different values.
		5. Repeat from 3 until there is no change in number of different vals
			Awesome!!

		"""

		# Some subroutines used in the canonicalization procedure

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
			"""Find the position of the maximum value in a list"""
			hi = li[0]
			hiPos = 0
			for i in range(len(li)):
				if li[i] > hi:
					hi = li[i]
					hiPos = i
				elif li[i] == hi:
					pass # TODO: Bond orders to differentiate?
			return hiPos

		matrix = self.connectMat # Shortcut XXX: Remove?

		size = len(matrix) # XXX: Remove this.

		# Atom labels -- used to build the new indices
		labels = [0 for x in range(size)]

		# Initially label each atom with its degree
		for i in range(size):
			for j in range(size):
				if matrix[i][j]:
					labels[i] += 1

		"""MORGAN'S ALGORITHM - PART ONE"""

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
			iterCount += 1 # XXX: Hack, "CCCCC" alkane gets stuck in inf loop

		"""MORGAN'S ALGORITHM - PART TWO"""

		# Score / Order the atoms by their labels.
		canonicalOrder = []

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

		"""RECONSTRUCT THE NEW MATRIX"""

		newMolMat = MolMatrix(self.size)

		def map_connectivity(newMat, newAtmPos, oldAtmPos, 
				newNbrPos, oldNbrPos):
			"""
			This function is called for each atom in order to map its
			previous old graph connectivity into the new one. (Atoms 
			are still bonded to the same neighbors.)
			"""
			def mapGraphs(newMat, oldMat):
				newMat[newAtmPos][newNbrPos] = \
						oldMat[oldAtmPos][oldNbrPos]
				newMat[newNbrPos][newAtmPos] = \
						oldMat[oldNbrPos][oldAtmPos]
			
			# Connectivity and Bond Order matrices
			mapGraphs(newMat.connectMat, self.connectMat)
			mapGraphs(newMat.bondOrderMat, self.bondOrderMat)


		def map_atom(newMat, newPos, oldPos):
			"""
			Transfer the atom type, charge, isotope, etc. data to its 
			new labeled position.
			"""
			newMat.atomTypes[newPos] = self.atomTypes[oldPos]
			newMat.atomCharges[newPos] = self.atomCharges[oldPos]
			newMat.atomIsotopes[newPos] = self.atomIsotopes[oldPos]

		# Build new matrix
		for i in range(size):
			newAtmPos = i
			oldAtmPos = canonicalOrder[i]
			neighbors = get_neighbors(matrix, oldAtmPos)

			map_atom(newMolMat, newAtmPos, oldAtmPos)

			for n in neighbors:
				oldNbrPos = n
				newNbrPos = canonicalOrder.index(n)

				map_connectivity(newMolMat, newAtmPos, oldAtmPos, 
							newNbrPos, oldNbrPos)

		return newMolMat

class ConstMolMatrix(MolMatrix):

	def __init__(self, mat):
		"""
		Build a constant molecular matrix from its predecessor.
		Once we've built the matrix, there is no reason we should
		modify it.
		"""

		def immutable(table):
			"""Convert a list-based matrix or connection table into an
			immutable tuple of tuples."""
			if not table:
				return ()
			table = table[:]
			for i in range(len(table)):
				if type(table[i]) == str:
					continue
				table[i] = tuple(table[i])
			return tuple(table)

		def make_immutable(func):
			"""Function decorator version."""
			def wrap(arg1, arg2=None):
				if arg2:
					return immutable(func(arg1, arg2))
				return immutable(func(arg1))
			return wrap

		@make_immutable
		def generate_alpha_table(mat):
			"""Compute the alpha (direct neighbor) connection table 
			upfront from the connection matrix."""
			sz = len(mat)
			table = [[] for x in range(sz)]
			for i in range(sz):
				for j in range(sz):
					if i == j:
						continue
					if mat[i][j]:
						table[i].append(j)
			return table

		@make_immutable
		def generate_beta_table(alpha):
			"""Compute the beta (neighbor of neighbor) connection table
			upfront from the alpha table."""
			sz = len(alpha)
			table = [[] for x in range(sz)]
			for i in range(sz):
				for j in alpha[i]:
					for k in alpha[j]:
						if k == i:
							continue
						table[i].append(k)
			return table

		def compute_hybridizations(bondOrderMat, neighborTable):
			"""
			Generate Hybridization State for each atom.
			This is determined by analyzing the number of pi bond
			systems from the connection table. 

			Each atom is one of: {'sp', 'sp2', 'sp3', 'error'}
			"""
			HYBRID_VALUES = {0: 'sp3', 1: 'sp2', 2: 'sp'}

			sz = len(bondOrderMat)
			hybrids = ['error' for x in range(sz)]
			for i in range(sz):
				numPi = 0 # Number of pi systems
				for n in neighborTable[i]:
					bond = bondOrderMat[i][n]
					if bond >= 2:
						numPi += bond - 1
				if numPi in HYBRID_VALUES:
					hybrids[i] = HYBRID_VALUES[numPi]

			return tuple(hybrids)

		def compute_degrees(atomTypes, neighborTable):
			"""
			Compute the atom degrees. For carbon, this is the number of
			other carbons it is directly attached to. For other atoms,
			it is the degree of the carbon it is attached to.
			"""
			# FIXME: Degree of ethers, amides, esters, carbonyls?
			# Granted these are 'functional groups' and not lone atoms. 
			def carbon_degree(atom):
				deg = 0
				for n in neighborTable[atom]:
					if atomTypes[n].upper() == 'C':
						deg += 1
				return deg

			sz = len(atomTypes)
			degrees = [-1 for x in range(sz)]
			for i in range(sz):
				aType = atomTypes[i].upper()
				deg = 0
				if aType == 'C':
					deg = carbon_degree(i)
				else: 
					# Calculate degree for non-carbon atoms
					deg = -1
					for n in neighborTable[i]:
						# FIXME: Not sure what to do if two carbons, 
						# eg. the Oxygen in 'COC'. For now, just use 
						# the first carbon.
						if atomTypes[n].upper() == 'C':
							deg = carbon_degree(n) # FIXME: Redundant
							break
				degrees[i] = deg
			return tuple(degrees)

		"""
		CLASS MEMBERS
			All matrices, connection tables, and label lists are are
			immutable tuples.
		"""

		# Number of atoms.
		self.size = mat.size

		# Connectivity matrix 
		self.connectMat = immutable(mat.connectMat)

		# Bond orders between atom pairs: 1, 1.5 (aromatic), 2, and 3
		self.bondOrderMat = immutable(mat.bondOrderMat)

		# Build neighbor/neighbor-of-neighbor tables. 
		self.alphaAtoms = generate_alpha_table(mat.connectMat)
		self.betaAtoms = generate_beta_table(self.alphaAtoms)

		# Atom labels (lists)
		self.types = tuple(mat.atomTypes) # C, O, N, Cl, etc.
		self.charges = tuple(mat.atomCharges)
		self.isotopes = tuple(mat.atomIsotopes)
		self.degrees = compute_degrees(mat.atomTypes, self.alphaAtoms)
		self.hybridizations = compute_hybridizations(mat.bondOrderMat, 
								self.alphaAtoms) # sp, sp2, sp3, or 'error'

	def print_matrix(self):
		"""
		Print the molecular data for debugging.
		"""
		print "Types: %s" % str(self.types)
		print "Charges: %s" % str(self.charges)
		print "Isotopes: %s" % str(self.isotopes)
		print "Hybridizations: %s" % str(self.hybridizations)
		print "Degrees: %s" % str(self.degrees)

		# XXX: Won't print >= 100 atoms nicely. Not that it would be
		# wise to print out such systems in the terminal...
		line1 = " "*5 if self.size < 10 else " "*6
		line2 = " "*5 if self.size < 10 else " "*6

		# Header atoms and header numbers
		for i in range(self.size):
			if len(self.types[i]) > 1:
				line1 += self.types[i]
			else:
				line1 += "%s " % self.types[i]
			if i < 10 or i %2 == 0:
				line2 += "%d " % i
			else:
				line2 += " "

		print line1
		print line2

		def row_header(i):
			# FIXME: Can't print systems with > 99 atoms or print 
			# 3 char atoms... but need for either is very unlikely
			atom = self.types[i] 
			if len(atom) < 2:
				atom += " "
			if self.size < 10 or i >= 10:
				return "%s %d  " % (atom, i)
			return "%s %d   " % (atom, i) 

		# Graph data
		for i in range(self.size):
			ln = row_header(i)
			for j in range(self.size):
				ln += str(int(self.bondOrderMat[i][j])) + " " \
						if self.bondOrderMat[i][j] else ". "
			print ln

		# Alpha atoms:
		print "\nAlpha Atoms:"
		for i in range(self.size):
			ln = row_header(i)
			ln += str(self.alphaAtoms[i])
			print ln

		# Beta atoms:
		print "\nBeta Atoms:"
		for i in range(self.size):
			ln = row_header(i)
			ln += str(self.betaAtoms[i])
			print ln


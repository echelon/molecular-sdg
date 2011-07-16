
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

 		# XXX: These values are cached.
		self._atomHybridizations = [None for x in range(size)]
		self._atomDegrees = [None for x in range(size)]
		self._neighbors = [None for x in range(size)]

	def numAtoms(self):
		"""Reports the number of atoms in the molecule."""
		# FIXME: Hydrogen exclusion, yet the ability to specify [nH] etc.
		# makes this inconsistent (and inaccurate--we can't get mol. weight)
		return self.size

	def getType(self, atomNum):
		return self.atomTypes[atomNum]

	def getNeighbors(self, atomNum):
		"""
		Returns the neighbors for a given atom. (ie. the alpha atoms).
		This data is gleaned from the connection matrix.
		"""
		if self._neighbors[atomNum] != None:
			return self._neighbors[atomNum][:]

		# TODO: Could cache this...
		n = []
		for i in range(self.size):
			if self.connectMat[atomNum][i]:
				n.append(i)

		self._neighbors[atomNum] = n[:]
		return n

	def getBetaAtoms(self, atomNum):
		"""
		Returns the beta atoms for a given atom (ie. the neigbors of
		neighbors). This data is gleaned from the connection matrix.
		"""
		# TODO: Could cache this...
		neighbors = self.getNeighbors(atomNum)
		beta = []
		for n in neighbors:
			for x in self.getNeighbors(n):
				if x != atomNum and x not in beta:
					beta.append(x)
		return beta

	def getHybridization(self, atomNum):
		"""
		Get the hybridization state of an atom.
		This is determined by analyzing the number of pi bond systems
		from the connection table. 

		Returns one of: {'sp', 'sp2', 'sp3', 'error'}
		"""
		# FIXME: Move to a higher level class?
		# FIXME: Generation and caching strategy could result in errors. 

		def generate():
			"""Generate Hybridization State for each atom in the 
			matrix."""
			hybridizations = {0: 'sp3', 1: 'sp2', 2: 'sp'}
			for i in range(self.numAtoms()):
				numPi = 0 # Number of pi systems
				for n in self.getNeighbors(i):
					bond = self.bondOrderMat[i][n]
					if bond >= 2:
						numPi += bond - 1
				hybrid = 'error'
				if numPi in hybridizations:
					hybrid = hybridizations[numPi]

				self._atomHybridizations[i] = hybrid

		if self._atomHybridizations[0] == None:
			generate()

		return self._atomHybridizations[atomNum]

	def getDegree(self, atomNum):
		"""
		Get the degree of an atom.
		For carbon, this is the number of other carbons it is directly
		attached to. For other atoms, it is the degree of the carbon it
		is attached to.
		"""
		# FIXME: Degree of ethers, amides, esters, carbonyls -- ???
		if self._atomDegrees[atomNum] != None:
			return self._atomDegrees[atomNum]

		aType = self.atomTypes[atomNum].upper()
		neighbors = self.getNeighbors(atomNum)
		deg = 0

		if aType == 'C':
			# For carbon atoms
			for n in neighbors:
				if self.atomTypes[n].upper() == 'C':
					deg += 1
		else: 
			# For non-carbon atoms
			deg = -1
			for n in neighbors:
				# FIXME: Not sure what to do if two carbons, eg. Ether. 
				if self.atomTypes[n].upper() == 'C':
					deg = self.getDegree(n)
					break

		self._atomDegrees[atomNum] = deg
		return deg

	def getConnectMat(self):
		"""Build a connection matrix we can use in Floyd's algorithm"""

		inf = float('Infinity')
		newMat = [[inf for x in range(self.size)] for xx in range(self.size)]

		for i in range(self.size):
			for j in range(self.size):
				val = self.connectMat[i][j]
				val = inf if not val else 1
				newMat[i][j] = val

		return newMat

	def print_matrix(self):
		"""Print the matrix. Debug."""
		print "TODO: Print Smiles Label, but from different module."

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


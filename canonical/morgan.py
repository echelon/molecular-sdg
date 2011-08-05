# TODO: No longer works. Needs to be refactored to support
# the current Molecule object.

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


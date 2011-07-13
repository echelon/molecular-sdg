#!/usr/bin/env python

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
		self.bondOrderMat = [[False for x in range(size)] 
									for xx in range(size)]

		# Atom Labels. 
		self.atomTypes = [False for x in range(size)]

	def print_matrix(self):
		"""Print the matrix. Debug."""
		print "TODO: Print Smiles Label, but from different module."
		#print "AdjMat for %s" % self.smiles
		# XXX: Won't print >= 100 atoms nicely. Not that I would want
		# to print out such systems in the terminal...
		ln = " "*3 if self.size < 10 else " "*4

		# Header atoms
		for i in range(self.size):
			if len(self.atomTypes[i]) > 1:
				ln += self.atomTypes[i]
			else:
				ln += "%s " % self.atomTypes[i]

		print ln
		ln = " "*3 if self.size < 10 else " "*4

		# Header numbers
		for i in range(len(self.connectMat)):
			if i < 10 or i %2 == 0:
				ln += "%d " % i
			else:
				ln += " "
		print ln

		# Graph data
		for i in range(len(self.connectMat)):
			if self.size < 10 or i >= 10:
				ln = "%d  " % i
			else:
				ln = "%d   " % i 
			for j in range(len(self.connectMat)):
				ln += str(int(self.connectMat[i][j])) + " " \
						if self.connectMat[i][j] else ". "
			print ln

	def numAtoms(self):
		"""Reports the number of atoms in the molecule."""
		return len(self.connectMat)

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

		# Build new matrix
		for i in range(size):
			newAtmPos = i
			oldAtmPos = canonicalOrder[i]
			neighbors = get_neighbors(matrix, oldAtmPos)

			newMolMat.atomTypes[newAtmPos] = self.atomTypes[oldAtmPos]

			for n in neighbors:
				oldNbrPos = n
				newNbrPos = canonicalOrder.index(n)

				map_connectivity(newMolMat, newAtmPos, oldAtmPos, 
							newNbrPos, oldNbrPos)

		return newMolMat


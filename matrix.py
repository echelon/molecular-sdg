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
		#print "AdjMat for %s" % self.smiles
		print "TODO: Label."
		# Won't print > 100 atoms nicely
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
				ln += str(self.connectMat[i][j]) + " " \
						if self.connectMat[i][j] else ". "
			print ln

	def numAtoms(self):
		"""Reports the number of atoms in the molecule."""
		return len(self.connectMat)

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

		matrix = self.connectMat

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
			iterCount += 1 # XXX: Hack, "CCCCC" alkane gets stuck in inf loop


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
		self.connectMat = newMatrix

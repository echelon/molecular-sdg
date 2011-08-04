"""
Ring Perception

TODO: BASIC DOCUMENTATION.

This is a somewhat outdated approach based on Zamora [1]. Despite
reading a review article [2], this is the only one I could find free
access to. I would really like to have [3] if anyone reading this would
like to donate a copy.

[1] Zamora, A., "An Algorithm for Finding the Smallest Set of Smallest
Rings", J. Chem. Inf. Comput. Sci., vol 16 (1976).

[2] Downs, G. M., et al., "Review of Ring Perception Algorithms for
Chemical Graphs", J. Chem. Inf. Comput. Sci., vol 29 (1989).

[3] Figueras, J. "Efficient Exact Solution of the Ring Perception
Problem", J. Chem. Inf. Comput. Sci., vol 34 (1994).
"""

# TODO: Implement phase 1 heuristics
# TODO: Implement phase 2
# TODO: Implement phase 3

def identify_rings(mol):
	"""
	Identify the rings in the system with Zamora's SSSR ring perception
	algorithm. 
	TODO: Further Documentation.
	"""

	def zamora_connectivity(mol):
		"""
		Zamora defines a connectivity index for each atom in a 
		molecular graph. This numerical index determines how crowded or
		connected each atom is, allowing us to start later processing
		with more central atoms. 
		"""
		# Mapping of number of neighbors to connectivity values
		# This is a huristic value assigned to each atom.
		# TODO: Valency > 4 is not handled by Zamora
		ki_values = {
			0: 0,
			1: 0,
			2: 1,
			3: 8,
			4: 64,
		}

		# Calculate each atom's connectivity
		ki = [0 for x in range(mol.size)]
		for atom in range(mol.size):
			n = len(mol.alphaAtoms[atom])
			ki[atom] = ki_values[n]

		# Now calculate each atom's neighborhood.
		# The congestedness of an atom's neighborhood yeilds higher
		# values for li. 
		li = [0 for x in range(mol.size)]
		for atom in range(mol.size):
			ksum = 0
			for neighbor in mol.alphaAtoms[atom]:
				ksum += ki[neighbor]

			li[atom] = ksum # Final heuristic score is the sum of ki and li.

		# However, ki is given an additional weight. 
		ci = [0 for x in range(mol.size)]
		for atom in range(mol.size):
			ci[atom] = 64* ki[atom] + li[atom]

		return ci

	def phase1(mol):
		"""
		Phase One processing takes care of "Type One" ring systems, as
		well as ring systems that can be reduced to Type One. 
		TODO: Documentation.
		"""

		# Connectivity index for each atom
		connectivity = zamora_connectivity(mol)

		# List of unused atom labels. (Elements deleted as used.)
		unusedAtoms = range(mol.size) # XXX: May not be best approach.

		rings = []
		while len(unusedAtoms) > 0:
			# Find the unused atom with the highest connectivity, then the 
			# smallest (and best) ring that contains it. 
			# TODO/FIXME/XXX: Absolutely need to implement Zamora's heuristics
			startAtom = connectivity.index(max(connectivity))
			ring = _find_smallest_ring(mol, startAtom)

			if len(ring) < 1:
				# Atom is not in a ring system. 
				unusedAtoms.remove(startAtom)
				connectivity[startAtom] = -1
				continue

			rings.append(ring)

			# Remove used atoms.
			unusedAtoms = filter(lambda x: x not in ring, unusedAtoms)

			# Negate connectivity for used atoms. 
			for i in range(len(connectivity)):
				if i not in unusedAtoms:
					connectivity[i] = -1

		return rings

	# TODO: Implement phase two and three.
	rings = phase1(mol)

	# Make constant tuples.
	for i in range(len(rings)):
		rings[i] = tuple(rings[i])

	return tuple(rings)

def _find_smallest_ring(mol, first=0, second=None):
	"""
	Find the smallest ring containing the atom 'first', or to find the
	smallest ring containing an edge, specify both 'first' and 'second'
	as atoms that are bonded together.

	Inputs:
		mol - molecule object w/ adj matrix, etc.
		first - The first atom of the path.
		second - (optional) The second atom of the path. In this case,
				 we desire to find the smallest ring containing the
				 edge (first, second).
	
	This is a dramatic rewrite of the ring finder algorithm in 
	Zamora [1]. As it was in PL/I, I couldn't tell that it was just DFS
	until I had collapsed all the evil GOTO statements.

	Ultimately, DFS could be replaced with another algorithm.
	"""

	class RingDFS(object):
		"""
		Contains metadata for a DFS search for rings.
		Interface aids in push()/pop() and check for search completion.
		"""

		def __init__(self, mol):
			"""CTOR"""
			sz = mol.size

			self.curAtom = -1	# Current Atom
			self.pathStack = []	# Current Path (DFS Stack)
			self.inPath = [False for x in range(sz)] # Fast per-atom lookup

			# Keep track of the children visited for each node. 
			# FIXME: A list of stacks, and thus inefficient to rebuild.
			self.neighbors = [[] for x in range(sz)]

			self._mol = mol

		def push(self, v):
			"""Push the atom onto the path."""
			self.pathStack.append(v)
			self.inPath[v] = True
			self.curAtom = v
			self.neighbors[v] = list(self._mol.alphaAtoms[v])

		def pop(self):
			"""Pop the atom off the path."""
			n = self.pathStack.pop()
			self.inPath[n] = False
			self.curAtom = self.pathStack[-1]
			self.neighbors[n] = list(self._mol.alphaAtoms[n])

		def foundCycle(self, n):
			"""If n is in the path, we found a cycle."""
			return self.inPath[n]

		def solutionFound(self):
			"""Whether the solution has been found."""
			# TODO: Verify this works.
			return len(self.pathStack) == 1 and \
					len(self.neighbors[self.curAtom]) == 0

	# ========== Begin Function ============

	dfs = RingDFS(mol)

	# Size of the smallest ring found; solution set.
	smallestRing = mol.size + 1
	resultSet = []

	# Initialize for ring finding at an atom
	dfs.push(first)

	# Initialize for ring finding at an edge
	if second != None:
		dfs.push(second)
		dfs.neighbors[first] = [] # Position at end.

	# Perform DFS
	while True:
		if len(dfs.neighbors[dfs.curAtom]) == 0:
			dfs.pop()
			if dfs.solutionFound():
				return resultSet
			continue

		n = dfs.neighbors[dfs.curAtom].pop(0)

		if len(dfs.pathStack) > 1 and n == dfs.pathStack[-2]:
			# Neighbor was just visited! It was the grandparent.
			# Don't go backwards, just ignore.
			continue

		if dfs.foundCycle(n):
			if n != first:
				# It's not root, so this won't work for us. Ignore.
				continue

			if len(dfs.pathStack) < smallestRing:
				# Save the result. 
				resultSet = dfs.pathStack[:]
				smallestRing = len(resultSet)
				dfs.pop()
				if dfs.solutionFound():
					return resultSet
				
			# XXX/TODO: We can insert huristics here.
			# Insert heuristic criteria to compare rings of equal size.
			# If our current ring is the same size of the current smallest,
			# we may wish for one to win out over the other based on some
			# heuristic. (Heteroatoms? Who knows. I'll have to see.)

			dfs.pop()
			if dfs.solutionFound():
				return resultSet

			continue

		else:
			dfs.push(n)
			if len(dfs.pathStack) > smallestRing:
				dfs.pop()
				if dfs.solutionFound():
					return resultSet

			continue


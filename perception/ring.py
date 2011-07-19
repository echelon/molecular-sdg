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

def find_smallest_ring(mol, first=0, second=None):
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
			sz = mol.numAtoms()

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
			self.neighbors[v] = self._mol.getNeighbors(v)

		def pop(self):
			"""Pop the atom off the path."""
			n = self.pathStack.pop()
			self.inPath[n] = False
			self.curAtom = self.pathStack[-1]
			self.neighbors[n] = self._mol.getNeighbors(n) # XXX: Correct?!

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
	smallestRing = mol.numAtoms() + 1
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


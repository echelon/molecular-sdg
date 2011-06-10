#!/usr/bin/env python

from smiles import Smiles

# XXX: Tis will be messy until I finish porting the algorithm
def test(string):
	"""Convert a SMILES string into a adjacency matrix representation."""

	# First, we must convert the input string into a proper queue

	ATOMS = ['H', 'C', 'N', 'P', 'O', 'S', 'F', 'Cl', 'Br', 'I']
	ATOMS_UPPER = ['H', 'C', 'N', 'P', 'O', 'S', 'F', 'CL', 'BR', 'I']

	BRANCHING = '()' # Branch Start & End
	BONDS = '=#'	 # Double and Triple Bonds
	CONNECTIVE = '%' # Connectivity beyond '9' -- TODO NOT YET HANDLED.

	smiles = Smiles(string)
	queue = smiles.tokens
	molgraph = MolGraph(smiles)
	graph = molgraph.graph

	# Convert to Graph

	root = None
	curNode = None
	nextBond = None
	connected = [] # Keep track of connectivity

	def backtrack_find_branchpoint(curNode):
		"""Return the pointer where the branch started."""
		while curNode:
			if curNode.branchStart:
				return curNode
			curNode = curNode.parent

		return False

	def backtrack_find_connectpoint(curNode, num):
		"""Return the pointer where the branch started."""
		while curNode:
			if curNode.connect == num:
				return curNode
			curNode = curNode.parent

		return False

	# Symbol type tests
	def is_atom(sym): return sym.upper() in ATOMS_UPPER
	def is_bond(sym): return sym in BONDS
	def is_branch_start(sym): return sym == '('
	def is_branch_end(sym): return sym == ')'
	def is_closure(sym): sym.isdigit()

	# Make note of the connection beween two atoms. 
	def connect(a1, a2):
		graph[a1][a2] = True
		graph[a2][a1] = True

	# Keep track of state during parsing
	atomCnt = 0 # Total num atoms
	atomPrev = 0 # Previous atom 
	branchStack = [] # Branching, eg. C(C)(C)C

	for sym in queue:

		if is_atom(sym):
			atomId = atomCnt
			atomCnt += 1

			if atomId != 0:
				connect(atomPrev, atomId)

			atomPrev = atomId

		elif is_branch_start(sym):
			branchStack.append(atomPrev)

		elif is_branch_end(sym): 
			atomPrev = branchStack.pop()


	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE
	# TODO: RESUME WORK HERE

	molgraph.print_graph()

	"""
	# TODO: 'cur' -> 'sym'
	# XXX: Abandoning this... start over above.
	atomCnt = 0
	prevAtomNum = 0
	for sym in queue:

		if sym.upper() in ATOMS_UPPER:
			# ATOMS
			# TODO: Record labels.

			if curNode:
				curNode.attachChild(newNode)

			curNode = newNode

			if nextBond:
				curNode.bond = nextBond
				nextBond = None

			prevAtomNum = atomCnt
			atomCnt += 1

		elif cur in BRANCHING:
			# BRANCHING
			if cur == '(':
				curNode.branchStart = True

			elif cur == ')':
				curNode = backtrack_find_branchpoint(curNode)
				curNode.branchStart = False

		elif cur in BONDS:
			nextBond = cur

		elif cur.isdigit():
			# CONNECTIVITY
			cur = int(cur)
			if cur not in connected:
				connected.append(cur)
				curNode.connect = cur
			else:
				oldNode = backtrack_find_connectpoint(curNode, cur)
				oldNode.connectsTo.append(curNode)
				curNode.connect = cur
				curNode.connectsTo.append(oldNode)

	"""


class MolGraph(object):
	"""Adjacency Matrix for molecules."""

	def __init__(self, molInput):
		"""CTOR"""

		if type(molInput) == str:
			molInput = Smiles(molInput)

		self.smiles = molInput
		
		sz = molInput.num_atoms()

		self.graph = [[False for x in range(sz)] for xx in range(sz)]


	def print_graph(self):
		"""Print the graph. Debug."""
		print "Graph for %s" % self.smiles
		# Won't print > 100 atoms nicely
		ln = " "*3 if self.smiles.num_atoms() < 10 else " "*4

		# Header numbers
		for i in range(len(self.graph)):
			if i < 10 or i %2 == 0:
				ln += "%d " % i
			else:
				ln += " "
		print ln

		# Graph data
		for i in range(len(self.graph)):
			if self.smiles.num_atoms() < 10 or i >= 10:
				ln = "%d  " % i
			else:
				ln = "%d   " % i 
			for j in range(len(self.graph)):
				ln += "1 " if self.graph[i][j] else ". "	
			print ln

# Temporary, for testing!!!
if __name__ == '__main__':
	import sys

	if len(sys.argv) < 2:
		print "Need to supply SMILES text as argument."
		sys.exit()
		
	mol = MolGraph(sys.argv[1])
	test(sys.argv[1])

	#mol.print_graph()	

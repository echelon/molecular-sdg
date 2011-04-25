import math

"""
Subset of Atoms used:
H
C
N P
O S
F Cl Br I
"""

class Atom(object):

	def __init__(self, ch):
		self.atom = ch
		self.children = []
		self.parent = None

		# State used to construct the graph's branching
		self.branchStart = False

		# Marker for cyclic systems
		self.connect = 0
		self.connectsTo = []

		# Bond type.
		self.bond = '-'
		self.isRing = False # TODO, also not in use.

		# Position Data
		self.x = 0
		self.y = 0

		# Orientation: 1 for up, -1 for down
		# TODO: This seems very primative. 
		self.orient = None # 1 or -1 

		# Empty spot availability
		self.available = {
			'left': True,
			'right': True,
			'up': True,
			'down': True,
		}

	def attachChild(self, node):
		node.parent = self
		self.children.append(node)

	def numAtoms(self):
		"""Returns the number of bonded atoms, MINUS hydrogens."""
		num = 0
		if(self.parent):
			num = 1
		num += len(self.children)
		return num

	def numChildren(self):
		return len(self.children)

	def geometry(self):
		"""Returns the geometry of the atom: linear, trigonal, tetrahedral."""
		bonds = []
		if self.parent:
			bonds.append(self.bond)
		for ch in self.children:
			bonds.append(ch.bond)

		double = bonds.count('=')
		triple = bonds.count('#')

		if triple == 1 or double == 2:
			return 'linear'
		elif double == 1:
			return 'trigonal'
		
		return 'tetrahedral'

	def hybridization(self):
		"""Returns a numerical value according to the hybridization of
		the atom: {0: 'lone atom', 1: 's', 2: 'sp', 3: 'sp2', 4: 'sp3'}
		"""
		count = len(self.children)
		if self.parent:
			count += 1

		if self.atom.upper() == 'C':
			# TODO: Check bond type...
			return 4

		return count

	def __str__(self):
		"""String Representation"""
		def make_string(atom, level=0):
			st = '\t' * level
			st += atom.bond + ' ' + atom.atom
			if atom.connect:
				st += '.' + str(atom.connect)
			st += '\n'
			for child in atom.children:
				st += make_string(child, level+1)
			return st
		return make_string(self)


	def positionLeft(self, mol):
		"""Position a molecule to the left of the current one."""
		angle = 180 - (109/2 + 90)
		mol.orient = self.orient * -1 # Invert orientation.

		endX = self.x + 30 # XXX: the length constant
		endY = self.y + math.tan(angle) * 4 * mol.orient * -1

		mol.x = endX
		mol.y = endY
	
		self.available['left'] = False
		mol.available['right'] = False

	def positionUp(self, mol):
		"""Position a molecule"""

		mol.orient = self.orient * -1 # Invert orientation.

		endX = self.x
		endY = self.y - 30 # XXX: Constant

		mol.x = endX
		mol.y = endY
	
		self.available['up'] = False
		mol.available['down'] = False

	def positionDown(self, mol):
		"""Position a molecule"""

		mol.orient = self.orient * -1 # Invert orientation.

		endX = self.x
		endY = self.y + 30 # XXX: Constant

		mol.x = endX
		mol.y = endY
	
		self.available['down'] = False
		mol.available['up'] = False

def create_graph(smileStr):

	# First, we must convert the input string into a proper queue

	ATOMS = ['H', 'C', 'N', 'P', 'O', 'S', 'F', 'Cl', 'Br', 'I']
	ATOMS_UPPER = ['H', 'C', 'N', 'P', 'O', 'S', 'F', 'CL', 'BR', 'I']

	BRANCHING = '()' # Branch Start & End
	BONDS = '=#'	 # Double and Triple Bonds
	CONNECTIVE = '%' # Connectivity beyond '9'

	READ_AHEAD = {
			'C': 'Cl', 
			'B': 'Br'
	}

	def make_queue(smileStr):
		"""Make a proper input queue from the raw string."""
		queue = []
		pos = 0
		while pos < len(smileStr):
			ch = smileStr[pos]

			# Check for Br, Cl, etc.
			if ch in READ_AHEAD and pos < len(smileStr)-1:
				rd = ch + smileStr[pos+1]
				if rd == READ_AHEAD[ch]:
					ch = rd
					pos += 1

			# TODO: Connectivity beyond '9'

			queue.append(ch)
			pos += 1

		return queue

	queue = make_queue(smileStr)

	# Convert to Minimum Spanning Tree

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

	for cur in queue:

		if cur.upper() in ATOMS_UPPER:
			# ATOMS
			newNode = Atom(cur)
			if not root:
				root = newNode

			if curNode:
				curNode.attachChild(newNode)

			curNode = newNode

			if nextBond:
				curNode.bond = nextBond
				nextBond = None

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

	return root 



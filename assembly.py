"""
Assembly Phase
TODO: Doc
"""

from math import * 
from chain import CHAIN_ANGLE
from ring import Point # TODO: Move elsewhere

def assemble(mol):
	"""
	TODO DOC
	This is the main exportable assembly process.
	"""
	queue = []
	headAtom = 0 # TODO: Heuristic selection (congestedness)

	queue.append(headAtom)
	isHeadAtom = True

	mol.pos[headAtom] = Point(0, 0)

	# Dequeue and draw seed. 
	while len(queue) > 0:
		seedAtom = queue.pop(0)
		substituent_placement(mol, seedAtom, isHeadAtom, queue)
		isHeadAtom = False


# XXX: Algorithm 6
# TODO: (1)  PFU Angular demand
# FIXME/VERIFY: (1) use of CFS[lo] vs hi
def angular_spacing(mol, seedAtom, isHeadAtom=False):
	"""
	Yeilds the optimum angular spacing between the remaining unplaced
	substituents. From [Helson]. 
	"""
	# With respect to the seed atom's substituents:

	angularDemand = 0 # Area taken up by substituent PFUs
	numSubstituent = 0 # Number of non-PFU substituents 

	# Bonds that have been placed 
	# Bond is 'unplaced' if either of its atoms is unplaced. 
	#placedBonds = [] 

	p = [] # Set of adjacent bonds (in PFUs) that have been encountered

	# Work with adjacent, unplaced bonds.
	for n in mol.alphaAtoms[seedAtom]:
		bond = frozenset([seedAtom, n])
		inPfu = False

		if mol.isPlaced[n] and mol.isPlaced[seedAtom]:
			continue

		if bond in p:
			continue

		# FIXME: Inefficient
		for ring in mol.rings:
			if bond in ring:
				# Increase by PFU demand, which is the complement of
				# CFS at atom n
				angularDemand += 0 # TODO TODO TODO TODO

				for bond in ring.bonds:
					p.append(bond)

				inPfu = True
				break
			
		if inPfu:
			numSubstituent += 1

	if isHeadAtom:
		numSubstituent -= 1

	# Correct for chain fixed angle spacing. 
	# Only if seed is a core chain atom (not end cap). 
	if mol.isInChain[seedAtom]:
		numSubstituent -= 1
		angularDemand += CHAIN_ANGLE

	# FIXME: Assume CFS lo.
	# Return final angular spacing calculation
	ang = 0
	try:
		ang = (mol.cfs[seedAtom]['lo'] - angularDemand) / (numSubstituent + 1)
	except:
		pass # TODO: Math exception handle

	return ang

# XXX: Algorithm 7
# TODO TODO TODO TODO: Entire algorithm.
def substituent_sequence(mol, seedAtom, isHeadAtom=False):

	return mol.alphaAtoms[seedAtom] # TODO TODO TODO TODO

	pass

	remaining = []
	complete = []

	# Return sequence
	sequence = []

	# TODO: Initialize to the set of atoms that have and have not been placed.
	remaining = list(mol.alphaAtoms[seedAtom])

	while len(remaining) > 0:

		s = None

		# TODO
		chain = False # XXX TEMP
		if seedAtom is chain:
			pass # TODO

		if len(complete) == 0:
			pass # TODO

		else:
			#s = substituent immediately CCW of last chosen atom
			pass # TODO

		"""
		ensure that s is not from the wrong side of a PFU --
			let b be the bond from the seed atom to s.
			if b is in a PFU:
				examine the PFU's bonds that are adj to the seed atom,
				and locate the bond q that is equal to the seed atom's local
				CFS.hi (there must be one).
				set s to q's other atom. 
				add the PFU's atoms to complete and subtract the PFU's 
				atoms from  remaining
		"""

		complete.append(s)
		remaining.remove(s)
		sequence.append(s)

	return sequence

# FIXME: Sequencing doesn't work yet...
# TODO: Ring case -- deposit PFU if bond is in PFU (& incr A.D.)
# TODO: Chain -- 'rightward bend' & double bond sterochemistry
# TODO: ENTIRE 'abnormal' case: Algo 10.
def substituent_placement(mol, seedAtom, isHeadAtom=False, drawQueue=[]):
	"""
	Place the substituents of a seed atom.
	Algorithms 8, 9, and 10 from [Helson].
	"""

	def ring_subst_placement():
		"""
		Substituent placement for ring seed atom
		Algorithm 8
		"""
		print "ring_subst_placement(%d)" % seedAtom

		beta = angular_spacing(mol, seedAtom, isHeadAtom)
		seq = substituent_sequence(mol, seedAtom, isHeadAtom)

		for n in seq:
			mol.cfs[seedAtom]['lo'] += beta

			place_atom(mol, seedAtom, n, False, drawQueue)

			# TODO: CHECK IF BOND IS IN PFU, IF SO DEPOSIT ENTIRE PFU
			# XXX XXX: Some conflict with algo11/placeAtom()
			#bond = frozenset([seedAtom, n])

	def chain_subst_placement():
		"""
		Substituent placement for core chain seed atom
		Algorithm 9
		"""
		print "chain_subst_placement(%d)" % seedAtom

		beta = angular_spacing(mol, seedAtom, isHeadAtom)
		seq = substituent_sequence(mol, seedAtom, isHeadAtom)

		chain = None
		for ch in mol.chains:
			if seedAtom in ch:
				chain = ch
				break

		for subst in seq:
			if mol.isPlaced[subst]:
				continue

			# TODO 'if takes rightward bend...'
			# If the substituent is in the same chain...
			rightward = True
			if subst in chain and rightward:
				mol.cfs[seedAtom]['lo'] += CHAIN_ANGLE
			else:
				mol.cfs[seedAtom]['lo'] += beta

			place_atom(mol, seedAtom, subst, False, drawQueue)

			# TODO: Double bond stereochemistry. 

	# XXX: Algorithm 10
	def uncat_subst_placement():
		#print "uncat_subst_placement(%d)" % seedAtom
		# TODO TODO TODO this is just a hack
		print "UNCATEGORIZED/TODO:"
		for n in mol.alphaAtoms[seedAtom]:
			place_atom(mol, seedAtom, n, False, drawQueue) # XXX HACK
		pass

	# Draw depending on seed atom type
	if mol.isInRing[seedAtom]:
		ring_subst_placement()

	elif mol.isInChain[seedAtom]:
		chain_subst_placement()

	else:
		uncat_subst_placement()


def place_atom(mol, seedAtom, placeAtom, seedIsHeadAtom=False, drawQueue=[]):
	"""
	Place Atom
	Algorithm 11 from [Helson].
	"""

	# XXX: added this check myself.
	if mol.isPlaced[placeAtom]:
		return

	mol.isPlaced[placeAtom] = True

	bond = frozenset([seedAtom, placeAtom])

	if seedIsHeadAtom and not mol.cfsInitialized[seedAtom]:

		mol.cfsInitialized[seedAtom] = True

		mol.cfs[seedAtom]['hi'] = 0
		mol.cfs[seedAtom]['lo'] = 0

		# FIXME: inefficient
		for ring in mol.rings:
			if bond in ring:
				mol.cfs[seedAtom]['hi'] = ring.angle # TODO: Doesn't work
				mol.cfs[seedAtom]['lo'] = ring.angle # TODO: DOesn't work
				break

	# TODO TODO TODO TODO
	# TODO TODO TODO TODO
	# TODO TODO TODO TODO
	# TODO TODO TODO TODO
	# PLACE A AT BONDLEN DISTANCE @ ANGLE CFS[lo]

	BOND_LEN = 100

	seedPt = mol.pos[seedAtom]
	cfs = mol.cfs[seedAtom]['lo']

	x = seedPt.x + BOND_LEN * sin(radians(cfs))
	y = seedPt.y + BOND_LEN * cos(radians(cfs))

	mol.pos[placeAtom] = Point(x, y)

	# FIXME: Degrees or radians? 
	# Set up the newly placed atom's CFS
	g = mol.cfs[seedAtom]['lo']
	dg = g - 180 # Points backwards

	mol.cfs[placeAtom]['hi'] = dg
	mol.cfs[placeAtom]['lo'] = dg
	mol.cfsInitialized[placeAtom] = True

	# END TODO TODO TODO TODO
	# END TODO TODO TODO TODO
	# END TODO TODO TODO TODO


	drawQueue.append(placeAtom) # TODO: ??? Why add to redraw queue?

	# TODO: IF BOND IN PFU, DEPOSIT WHOLE PFU NOW.



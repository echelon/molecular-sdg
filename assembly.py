"""
Assembly Phase
TODO: Doc
"""

from math import * 
from chain import CHAIN_ANGLE
from point import Point
from util.angles import normalize

def assemble(mol):
	"""
	TODO DOC
	This is the main exportable assembly process.
	"""
	queue = []
	headAtom = 0 # TODO: Heuristic selection (congestedness)

	queue.append(headAtom)
	isHeadAtom = True

	print "Atom 0 is at (0, 0)"

	mol.pos[headAtom] = Point(0, 0)
	mol.isPlaced[headAtom] = True

	# Dequeue and draw seed. 
	while len(queue) > 0:
		seedAtom = queue.pop(0)
		substituent_placement(mol, seedAtom, isHeadAtom, queue)
		isHeadAtom = False


# XXX: Algorithm 6
# XXX XXX XXX XXX XXX: Something is VERY WRONG with this algorithm
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

	p = [] # Set of adjacent bonds (in PFUs) that have been encountered

	# Work with adjacent, unplaced bonds.
	for n in mol.alphaAtoms[seedAtom]:
		bond = frozenset([seedAtom, n])
		inPfu = False

		# Only 'placed' if both atoms are placed.
		if mol.isPlaced[n] and mol.isPlaced[seedAtom]:
			continue

		if bond in p:
			continue

		# Add angle demand of ring substituents.
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
			
		if not inPfu:
			numSubstituent += 1

	if isHeadAtom:
		numSubstituent -= 1

	# Correct for chain fixed angle spacing. 
	# Only if seed is a core chain atom (not end cap). 
	if mol.isInChain[seedAtom]:
		numSubstituent -= 1
		angularDemand += CHAIN_ANGLE

	ang = 0
	cfs = mol.cfs[seedAtom]['lo'] # FIXME: Assume CFS lo.
	try:
		ang = (cfs - angularDemand) / (numSubstituent + 1)
		return normalize(ang)
	except:
		print "Angular Spacing exception occurred... using second algo"
		pass # TODO: Math exception handle

	return 0


def angular_spacing2(mol, seedAtom, isHeadAtom=False):
	"""
	Attempts to fix problems with [Helson] algo.
	"""
	try:
		return normalize(360/len(mol.alphaAtoms[seedAtom]))
	except:
		return 360

# XXX: Algorithm 7
# TODO TODO TODO TODO: Entire algorithm.
def substituent_sequence(mol, seedAtom, isHeadAtom=False):
	"""
	Sequence the substituents into attachment order.
	"""
	alpha = mol.alphaAtoms[seedAtom]
	left = [x for x in alpha if mol.isPlaced[x] == False]
	done = [x for x in alpha if mol.isPlaced[x] == True]
	seq = []

	chain = None
	if mol.isInChain[seedAtom]:
		for ch in mol.chains:
			if seedAtom in ch:
				chain = ch
				break

	# Add chain members first (not in any kind of order...)
	# XXX XXX: This is NOT the algorithm from [Helson]
	if chain:
		i = chain.index(seedAtom)
		l = chain[i-1] if i > 0 else chain.caps[0]
		r = chain[i+1] if i < len(chain)-1 else chain.caps[1]
		if l in left:
			left.pop(left.index(l))
			done.append(l)
			seq.append(l)

		if r in left:
			left.pop(left.index(r))
			done.append(r)

	# XXX: Add remaining atoms. NOT BASED ON [Helson]!
	# TODO: By priority
	done.extend(left)
	seq.extend(left)

	# TODO: Which to return?
	return done

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
		print "\nring_subst_placement(%d)" % seedAtom

		beta = angular_spacing(mol, seedAtom, isHeadAtom)
		seq = substituent_sequence(mol, seedAtom, isHeadAtom)

		for n in seq:
			mol.cfs[seedAtom]['lo'] += beta
			mol.cfs[seedAtom]['lo'] = normalize(mol.cfs[seedAtom]['lo'])

			place_atom(mol, seedAtom, n, False, drawQueue)

			# TODO: CHECK IF BOND IS IN PFU, IF SO DEPOSIT ENTIRE PFU
			# XXX XXX: Some conflict with algo11/placeAtom()
			#bond = frozenset([seedAtom, n])

	def chain_subst_placement():
		"""
		Substituent placement for core chain seed atom
		Algorithm 9
		"""
		print "\nchain_subst_placement(%d)" % seedAtom

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

			cfslo = mol.cfs[seedAtom]['lo']
			print "\n\t* cfs.lo (before change): %f" % cfslo
			print "\t* beta angle: %f" % beta
			print "\t* FxAS chain angle: %f" % CHAIN_ANGLE

			# Handle chain zigzag.
			# TODO 'if takes rightward bend...'
			# If the substituent is in the same chain...
			rightward = False
			if chain and subst in chain:
				z = chain.index(seedAtom)
				if z < len(chain.zigzag) and chain.zigzag[z].upper() == 'R':
					rightward = True

			if chain and subst in chain:
				if rightward:
					print "\t* R zigzag"
					#mol.cfs[seedAtom]['lo'] += beta # XXX TEST
					mol.cfs[seedAtom]['lo'] += CHAIN_ANGLE
					mol.cfs[seedAtom]['lo'] = normalize(mol.cfs[seedAtom]['lo'])
				else:
					print "\t* L zigzag"
					#mol.cfs[seedAtom]['lo'] -= beta # XXX TEST
					mol.cfs[seedAtom]['lo'] -= CHAIN_ANGLE # XXX
					mol.cfs[seedAtom]['lo'] = normalize(mol.cfs[seedAtom]['lo'])

			else:
				print "\t* NOT IN CHAIN!"
				mol.cfs[seedAtom]['lo'] += beta # XXX TEMP
				mol.cfs[seedAtom]['lo'] = normalize(mol.cfs[seedAtom]['lo'])



			print "\t** NOW CFS.LO IS: %f" % mol.cfs[seedAtom]['lo']
			place_atom(mol, seedAtom, subst, False, drawQueue)

			# TODO: Double bond stereochemistry. 

	# XXX: Algorithm 10
	def uncat_subst_placement():

		print "\nuncat_subst_placement(%d)" % seedAtom

		# XXX XXX XXX : Custom algo -- will it fix the problems?
		angle = angular_spacing2(mol, seedAtom, isHeadAtom)
		seq = substituent_sequence(mol, seedAtom, isHeadAtom)

		# May be able to make pseudo trigonal planar!
		# TODO TODO TODO TODO Must also ensure no triple bonds and no
		# more than two single bonds. Need data structure in Molecule
		# to keep track of this. I don't want to parse here. 
		pseudoTrigonal = False
		if len(mol.alphaAtoms[seedAtom]) == 2: # TODO
			print "\t* Pseudotrigonal"
			pseudoTrigonal = True
			angle = 120
		else:
			print "\t* NOT Pseudotrigonal"

		for subst in seq:
			if mol.isPlaced[subst]:
				continue

			print "\n\t* AngleSpacing: %0.0f" % angle
			cfslo = mol.cfs[seedAtom]['lo']
			print "\t* cfs.lo = %f" % cfslo

			# Update CFS
			if not isHeadAtom:

				if pseudoTrigonal:
					# TODO: ADD OR SUBTRACT ANGLE DEPENDING ON ZIGZAG SENSE
					# TODO: some other case for Zigzag??
					pass
				else:
					print angle
					mol.cfs[seedAtom]['lo'] += angle
					mol.cfs[seedAtom]['lo'] = normalize(mol.cfs[seedAtom]['lo'])

			cfslo = mol.cfs[seedAtom]['lo']
			print "\t* cfs.lo = %f" % cfslo

			print ""
			place_atom(mol, seedAtom, subst, False, drawQueue)

		# TODO: Double check dbl bond stereochem

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

	# FIXME 1: Need to read bond length from global configs
	# FIXME 2: Switch to radians

	BOND_LEN = 100

	seedPos = mol.pos[seedAtom]
	cfs = mol.cfs[seedAtom]['lo']

	x = seedPos.x + BOND_LEN * cos(radians(cfs))
	y = seedPos.y + BOND_LEN * sin(radians(cfs))
	mol.pos[placeAtom] = Point(x, y)

	# Set up the newly placed atom's CFS (points backwards!)
	g = mol.cfs[seedAtom]['lo'] - 180
	mol.cfs[placeAtom]['hi'] = g
	mol.cfs[placeAtom]['lo'] = g
	mol.cfsInitialized[placeAtom] = True

	# XXX XXX DEBUG DUBUG
	pt = mol.pos[seedAtom]
	cfs1 = mol.cfs[seedAtom]['lo']
	cfs2 = mol.cfs[placeAtom]['lo']
	print "\tAtom %d Placed at %f:\n\t\t%s" % (placeAtom, cfs1, pt)
	print "\t\tseed CFS.lo = %f, \n\t\tplaced CFS.lo = %f" % (cfs1, cfs2)

	drawQueue.append(placeAtom) # TODO: ??? Why add to redraw queue?

	# TODO: IF BOND IN PFU, DEPOSIT WHOLE PFU NOW.
	# TODO: DEPOSIT RING GROUPS -- MUST CALCULATE OFFSET AND ANGLES
	# XXX: Not based on Helson!
	if mol.isInRing[placeAtom]:
		# Calculate the atom offset position in the ring group
		# TODO/FIXME: Doesn't take angle into consideration...
		offset = Point(0, 0)
		ringGroup = mol.ringGroupRef[placeAtom]
		for ring in ringGroup:
			if placeAtom in ring:
				absPos = mol.pos[placeAtom]
				relPos = ring.pos[ring.index(placeAtom)]
				offset.x = absPos.x - relPos.x
				offset.y = absPos.y - relPos.y
				break

		# Assign positions
		for ring in ringGroup:
			for atom in ring:
				pos = ring.pos[ring.index(atom)]
				pos.x += offset.x
				pos.y += offset.y
				mol.pos[atom] = pos
				mol.isPlaced[atom] = True

		# Queue all ring atoms...
		for ring in ringGroup:
			for atom in ring:
				drawQueue.append(atom)




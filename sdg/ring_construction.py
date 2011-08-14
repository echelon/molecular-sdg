"""
Ring Construction picks up where Ring Analysis left off.
It is the process of laying out the ring coordinates.
"""

from math import *
import sys

from ring import RING_TYPES, Point

# TODO: Bridged ring positioning
# TODO: Spiro ring positioning
# FIXME: cw/ccw direction handling (how does it work?)

def construct_group(ringGroup):
	"""
	Construct the coordinates, CFS, etc. for a ring group.
	Follows from [Helson] p~335
	"""
	remRings = list(ringGroup.peelOrder[:])

	# Take care of core ring
	core = remRings.pop()
	if core.type != RING_TYPES.CORE:
		raise Exception, "Last Ring from Ring Peeling is not Core!"

	# Assign points for core (using CW)
	positions = regular_polygon(core, bondLen=50.0, direc='cw')
	for i in range(len(positions)):
		core.pos[i] = positions[i]

	# Work on remaining rings.
	# TODO: Function, "attach_fused"

	assigned = [core]

	while remRings:
		ring = remRings.pop()

		# Get the fusion atoms
		fusionAtoms = None
		fusionRing = None
		for r in assigned:
			atoms = ring.bonds & r.bonds
			if not atoms:
				continue
			fusionAtoms = list(tuple(atoms)[0])
			fusionRing = r
			break

		a = fusionAtoms[0]
		b = fusionAtoms[1]

		# Copy already known fused-edge positions to the ring
		# They'll be used in the regular_polygon() procedure. 
		ring.pos[ring.index(a)] = fusionRing.pos[fusionRing.index(a)]
		ring.pos[ring.index(b)] = fusionRing.pos[fusionRing.index(b)]

		# XXX/TEMPORARY -- THIS IS JUST FOR DEBUG
		ring.fusionAtom[ring.index(a)] = True
		ring.fusionAtom[ring.index(b)] = True
		fusionRing.fusionAtom[fusionRing.index(a)] = True
		fusionRing.fusionAtom[fusionRing.index(b)] = True
		# END XXX/TEMPORARY -- THIS IS JUST FOR DEBUG

		# Switch atoms depending on which side is being drawn.
		# FIXME/TODO: Does not support 'ccw' drawing.
		#seq = ring.sequence('cw')
		#aPos = seq.index(a)
		#bPos = seq.index(b)
		#if (aPos + 1) % len(seq) != bPos:
		if ring.getDirection() == 'ccw':
			ring.swapped = True # XXX/TEMPORARY: BOOL IS DEBUG ONLY
			print "\t>> SWITCH DIRECTION"
			a, b = b, a

		positions = regular_polygon(ring, atomA=a, atomB=b, 
									bondLen=50.0, direc='cw')

		# Map positions.
		aPos = ring.index(a)
		bPos = ring.index(b)

		aPtPos = None
		bPtPos = None
		for i in range(len(positions)):
			if positions[i] == ring.pos[aPos]:
				aPtPos = i
			if positions[i] == ring.pos[bPos]:
				bPtPos = i

		try: 

			rDirec = False # Left
			if (aPos + 1) % len(ring) != bPos:
				rDirec = True # Right

			pDirec = False # Left
			if (aPtPos + 1) % len(positions) != bPtPos:
				pDirec = False # Right

		except:
			print "========================================"
			print "======== Mapping Positions to Ring ====="
			print "Fusion Atoms: %d and %d" % (a, b)
			print "Indices: %d and %d" % (aPos, bPos)
			print "Map Points: "
			print "\t", 
			print aPtPos
			print "\t", 
			print bPtPos
			print ""
			print "Positions Generated: "
			print positions
			print ""
			print "Supposed to have included:"
			print "Position A: %s" % ring.pos[aPos]
			print "Position B: %s" % ring.pos[bPos]
			print ""

			sys.exit()
			
		if rDirec == pDirec:
			for i in range(len(ring)):
				ri = (aPos + i) % len(ring)
				pi = (aPtPos + i) % len(ring)
				ring.pos[ri] = positions[pi]

		else:
			for i in range(len(ring)):
				ri = (aPos + i) % len(ring)
				pi = (aPtPos - i) % len(ring) # XXX: Note sign
				ring.pos[ri] = positions[pi]
		
		assigned.append(ring)

		#sys.exit()

	return

# TODO for regular_polygon: handle cw/ccw. Still not sure how it works.
def regular_polygon(ring, atomA=None, atomB=None, bondLen=100.0, direc='cw'):
	"""
	Calculate the coordinate positions for a regular polygon.

	Inputs:
		ring - ring object
		atomA/atomB - desired first and second atom
		bondLen - desired length of the bonds.
		direc - direction of the construction.

	Objectives: 
		1. Find center point O
		2. Find angle to first vertex
		3. Find vertex points by adding or subtracting phi 
		   once per vertex. (Easy).
	
	Returns: A list of Point object coordinates for each vertex.
	"""
	if direc not in ('cw', 'ccw'):
		raise Exception, 'Direction supplied `%s` is invalid.' % direc

	# TODO: Update to ensure proper cw/ccw handling

	ptA = None # Point objects
	ptB = None
	idxA = 0 # Ring index of the atom, that is ring[idx]
	idxB = 0

	# If atoms are not specified, this is probably a core ring, and we
	# are free to align the ring with the coordinate system. All other
	# rings should be looked up in the ring index.

	if atomA != None and atomB != None:
		idxA = ring.index(atomA)
		idxB = ring.index(atomB)
		ptA = ring.pos[idxA]
		ptB = ring.pos[idxB]
		"""
		print "=======================================\n"
		print "Atoms: %d, %d" % (atomA, atomB)
		print "Indices: %d, %d" % (idxA, idxB)
		print "Position 1: %s" % str(ptA)
		print "Position 2: %s" % str(ptB)
		print ""
		"""
	else:
		idxA = 0
		# TODO: is this correct cw/ccw handling?
		idxB = ((idxA + 1) if direc == 'cw' else (idxA - 1)) % len(ring)

		# Align to coordinate system.
		ptA = Point(0.0, 0.0)
		if len(ring) % 2 == 0:
			ptB = Point(ptA.x, ptA.y + bondLen)
		else: 
			ptB = Point(ptA.x + bondLen, ptA.y)

	# Discern bond length.
	# FIXME: Inconsistent to recalculate.
	x = abs(ptA.x - ptB.x)
	y = abs(ptA.y - ptB.y)
	L = sqrt(x**2 + y**2) # Bond length

	# Characteristic angle; the angle between two vertices (from
	# polygon center 'O'), eg. angle <AOB
	phi = radians(360.0)/len(ring)

	# Distance from polygon center 'O' to each vertex, eg. O->A
	r = L/(2*sin(phi/2))

	# Point between the line joining A and B
	cX = ptA.x + 0.5*(ptB.x - ptA.x)
	cY = ptA.y + 0.5*(ptB.y - ptA.y)
	ptC = Point(cX, cY)

	# Distance from polygon center 'O' to the point in the 
	# center of line AB, 'C'. ie. O->C
	z = sqrt(r**2 - 0.25*L**2) # Also: z = r*cos(phi/2)

	# Scaled Perpendicular direction vector from C to O
	# http://answers.google.com/answers/threadview/id/419874.html
	k = L/2.0
	u = ((ptA.y - ptC.y)/k,
		(ptC.x - ptA.x)/k)

	# Polygon central point, calculated with the direction vector.
	oX = ptC.x + z* u[0]
	oY = ptC.y + z* u[1]
	ptO = Point(oX, oY)

	# Find angle to point A. 
	adj = ptA.x - ptO.x
	opp = ptA.y - ptO.y

	theta = 0
	if adj == 0.0:
		# TODO: Verify
		print "Must use predefined angle." # XXX DEBUG
		if opp >= 0:
			theta = pi * 1/2
		else:
			theta = pi * 3/2
	else:
		theta = atan(opp/adj)

	if not atomA or not atomB: # XXX: DEBUG
		atomA = 0 # XXX DEBUG
		atomB = 0 # XXX DEBUG

	"""
	print "-----------"
	print "Old:"
	print "%s idx: %d atom: %d" % (str(ptA), idxA, atomA)
	print "%s idx: %d atom: %d" % (str(ptB), idxB, atomB)
	print "-----------"
	print ""
	print "Point O: %s" % str(ptO) # CORRECT
	print "R: %f" % r
	print ""
	"""

	# Calculate the atom positions in the ring.
	i = idxA
	positions = []
	for x in range(len(ring)):
		px = ptO.x + cos(theta) * r
		py = ptO.y + sin(theta) * r
		positions.append(Point(px, py))

		"""	
		print "Angle: %f" % degrees(theta)
		print "Point: %s" % str(ring.pos[i])
		print "Index: %d Atom: %d" % (i, ring[i])
		print ""
		"""

		i = (i + 1) % len(ring)
		theta += phi # XXX: CW

	return positions


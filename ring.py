"""
Ring-related data structures and essential functions.

	Classes:
		Ring
		RingGroup
		RingSystem (TODO: Not necessarily needed...)

	Functions:
		partition_rings()
"""

from util.enum import enum

"""
Ring types are discovered during ring peeling. As each ring is peeled
from the ring group, it is assigned one of the types below.
"""
RING_TYPES = enum(
	'CORE', 		# Last remaining ring (tractable case)
	'TOUGH_CORE',	# Remaining ring(s) (intractable case)
	'SPIRO',		# Spiro rings share one atom
	'FUSED', 		# Fused rings share an edge
	'BRIDGED', 		# Bridged rings share two non-adjacent atoms
	'IRREGULAR',	# Irregular ring
	'NONE'			# No type assigned
)

class Point(object):
	"""Represents a 2D position."""
	def __init__(self, x=None, y=None):
		self.x = x
		self.y = y
	def __str__(self):
		x = 'N' if self.x == None else str(self.x)
		y = 'N' if self.y == None else str(self.y)
		return "pt(%s, %s)" % (x, y)
	def __repr__(self):
		return str(self)

# FIXME: Remove most of Ring's methods. They belong in analysis only!
class Ring(tuple):
	"""
	Ring Data Structure.
	Contains the atom cycle in the ring, bonds, etc.
	Also assists in processing tasks.
	"""

	def __new__(cls, path, positions=None):
		"""
		Build tuple subclass instance. The underlying tuple holds the
		ring path itself, so we do some cleanup to ensure each atom
		only occurs once.

		The start and end atom are bonded. (Technically, rings are also
		undirected cycles. There is no real 'start' or 'end', or even
		'direction'.
		"""
		def build_path(path):
			"""Remove the end atom if it is also the start atom."""
			newPath = list(path[:])
			if newPath[0] == newPath[-1]:
				newPath.pop()
			return tuple(newPath)

		# TODO: Throw exception on twice-included atoms?
		return tuple.__new__(cls, build_path(path))

	def __init__(self, path, positions=None):
		"""
		Ring Constructor
		Must specify the atom cycle that constitutes the ring.
		"""
		def build_bonds(path):
			"""Build a tuple of sets that represent the bonds."""
			l = len(path)
			bonds = []
			for i in range(l):
				a = path[i]
				a2 = path[(i+1)%l]
				bonds.append(frozenset([a, a2]))
			return frozenset(bonds)

		# XXX/FIXME/TODO: DEPRECATED! This is handled by the tuple 
		# superclass.
		self.ringPath = self[:]

		# Bonds in the ring. 
		# A set of sets (each inner set is an edge), which makes testing
		# for bridged rings easy. 
		self.bonds = build_bonds(self.ringPath)

		# Assigned ring group and in-group ID/offset assigned by
		# RingGroup.
		self.group = None
		self.groupOffset = None

		# The type of ring (elucidated during ring peeling).
		# TODO: rename peelStrategy.
		self.type = RING_TYPES.NONE

		# Ring-Local coordinate position of every atom in the ring.
		self.pos = None

		if positions and len(self) == len(positions):
			self.pos = positions
		else:
			self.pos = [Point() for x in range(len(self[:]))]

		# Ring-Local CFS (in radians) for every atom in the ring.
		self.cfs = [{'hi':0, 'lo':0} for x in range(len(self[:]))]

	def getDirection(self):
		"""Determine which direction the ring is directed."""
		# FIXME: Somewhat expensive?
		a = self[0]
		b = self[1]

		seq = self.sequence('cw')
		ap2 = seq.index(a)
		bp2 = seq.index(b)

		if (ap2 + 1) % len(self) == bp2:
			return 'cw'
		return 'ccw'

	def sequence(self, direc='cw'):
		"""
		Reorder the ring clockwise (cw) or counter-clockwise (ccw).

		Input: direction ('cw' or 'ccw')
		Output: new Ring with:
			* reordered atoms
			* reordered Positions (references)
		"""
		# FIXME: Simplify
		if direc not in ['cw', 'ccw']:
			raise Exception, "Invalid sequence direction, `%s`." % direc

		atoms = self[:]
		sz = len(self)

		# Get atom with lowest Y position.
		lowY = None
		lowYPos = 0 # Defaults to 0
		for i in range(len(self)):
			y = self.pos[i].y
			if not y:
				continue
			if not lowY or y < lowY:
				lowY = y
				lowYPos = i

		i = lowYPos
		initialI = i
		path = []
		positions = []

		if direc == 'cw':
			while len(path) < sz:
				path.append(atoms[i])
				positions.append(self.pos[i])
				i = (i + 1) % sz # Next: Left atom
		else:
			while len(path) < sz:
				path.append(atoms[i])
				positions.append(self.pos[i])
				i = (i - 1) % sz # Next: Right atom

		return Ring(path, positions)

	# TODO: Deprecated. Move into relevant module.
	def isCentralRing(self, ringList):
		"""
		Determine ring centrality.
		If a ring is central, its removal will partition the remining
		rings. We want to peel the rings on the molecule extremety
		first.

		This algorithm is modified from [Helson] in order to compensate
		for molecules such as PhPh-PhPhPh.
		"""
		# TODO: Update to work with RingGroups instead. It'll be much faster.
		rings = list(ringList[:])
		numRings = len(rings)

		# Remove our ring from the system
		pos = None
		for i in range(len(rings)):
			if rings[i] == self:
				pos = i
				break

		if pos != None:
			rings.pop(pos)
			numRings -= 1

		# Arbitrarily choose starting ring. To compensate for an issue
		# with Helson's algoritm, we'll choose the first that is spiro
		# or fused to the current ring. 
		rootPos = None
		for i in range(len(rings)):
			if self.isSpiroTo(rings[i]) or self.isFusedTo(rings[i]):
				rootPos = i
				break

		if rootPos == None:
			rootPos = 0

		# Iteratively collect rings that are fused/spiro.
		# If the number collected equals the number of rings (not
		# including this one), then this ring is not "central". 
		collected = [rootPos]
		justAdded = [rootPos]
		while len(justAdded) != 0:
			remRings = [x for x in range(len(rings)) if x not in collected]
			added = []
			for r in justAdded:
				for s in remRings:
					if rings[r].isSpiroTo(rings[s]) or \
							rings[r].isFusedTo(rings[s]):
						collected.append(s)
						added.append(s)
			justAdded = added

		if len(collected) == numRings:
			return False

		# If not all rings were gathered, perhaps the ring system isn't
		# "connected". Perform the collection again, this time with 
		# this ring included. Compare the results. If the number of 
		# collected rings (minus this one) equals the same number as
		# before, then this ring had no "connecting" effect and is not
		# central.
		rings.append(self)

		collected2 = [rootPos]
		while len(justAdded) != 0:
			remRings = [x for x in rings if x not in collected]
			added = []
			for r in justAdded:
				for s in remRings:
					if r.isSpiroTo(rings[s]) or r.isFusedTo(rings[s]):
						collected.append(s)
						added.append(s)
			justAdded = added

		if len(collected2) - 1 == len(collected):
			return False

		return True

	# TODO: Deprecated. Move into relevant module.
	def isSpiroTo(self, otherRing):
		"""
		Returns True if the ring shares only one atom with the other
		ring.
		"""
		count = 0
		for a in self.ringPath:
			if a in otherRing.ringPath:
				count += 1
		return count == 1

	# TODO: Deprecated. Move into relevant module.
	def isFusedTo(self, otherRing):
		"""
		Returns True if the ring shares one bond (edge) with the other
		ring.
		"""
		count = 0
		for b in self.bonds:
			if b in otherRing.bonds:
				count += 1
		# TODO: Can they share more than one and be considered fused?
		return count == 1

	# TODO: Deprecated. Move into relevant module.
	def isBridgedTo(self, otherRing):
		"""
		Returns True if the ring bridges the other ring.
		"""
		# TODO: Bridged ring support
		#print "TODO: Bridged ring support"
		return False

	def __repr__(self):
		"""Returns representation of object."""
		# Type => String map. Don't handle 'none'
		rtype = ''
		if self.type == RING_TYPES.CORE:
			rtype = 'Core'
		elif self.type == RING_TYPES.TOUGH_CORE:
			rtype = 'ToughCore'
		elif self.type == RING_TYPES.SPIRO:
			rtype = 'Spiro'
		elif self.type == RING_TYPES.FUSED:
			rtype = 'Fused'
		elif self.type == RING_TYPES.BRIDGED:
			rtype = 'Bridged'
		elif self.type == RING_TYPES.IRREGULAR:
			rtype = 'Irregular'

		if rtype:
			return "%sRing%s"  % (rtype, str(self[:]))
		return "Ring%s" % str(self[:])

	def __str__(self):
		"""Returns string representation of the object."""
		st = repr(self)
		st += "\n Positions:\n"
		for i in range(len(self.pos)):
			st += "   * %d : %s\n" % (self[i], str(self.pos[i]))
		return st

class RingGroup(tuple):
	"""
	Container for all rings that are connected together.
	Eg. "PhPh-C-Ph-C-PhPhPh" contains 3 ring groups.

	See partition_rings() in this module for the function that creates
	ring groups.

	This is a tuple subclass, so Pythonic usage should work well:

	>>> group = RingGroup([ring1, ring2, ring3], rgId=1)
	>>> ring2 in group
	True
	>>> group[:]
	(ring1, ring2, ring3)

	RingGroups are useful for saving on computational legwork later: we
	can restrict searches or traversals to be within the group instead
	of searching through the entire molecule's ring system.
	"""

	def __new__(cls, rings, rgId = None):
		"""Build tuple subclass instance."""
		# TODO: Check to ensure no two rings are included twice.
		return tuple.__new__(cls, rings)

	def __init__(self, rings, rgId=None):
		"""
		Input:
		* a list of Ring objects (handled by __new__), and is the
		  tuple data
		* an identifier (optional)
		"""

		# An arbitrary ID that the grouping algorithm assigns. Should
		# be unique for every ring group in the molecule, but it is 
		# not mandatory.
		self.rgId = rgId

		# Peel order established in ring analysis.
		self.peelOrder = []

		# Connection tables. The subscripts and values used are NOT the
		# atom labels, rather they are the subscripts internal to the
		# ring.
		self.spiroTo = [[] for x in range(len(self))]
		self.fusedTo = [[] for x in range(len(self))]
		self.bridgedTo = [[] for x in range(len(self))]

		# Assign in-group offsets to each ring.
		for i in range(len(self)):
			rings[i].group = self
			rings[i].groupOffset = i

		# Build fused table.
		for i in range(len(self)):
			for j in range(len(self)):
				if i == j:
					continue
				if self[i].bonds & self[j].bonds:
					self.fusedTo[i].append(j)

	def __repr__(self):
		"""Debug representation."""
		if self.rgId == None:
			return "<RingGroup(%d rings)>" % len(self)
		return "<RingGroup(#%d, %d rings)>" % (self.rgId, len(self))

	def __str__(self):
		"""String representation."""
		ret = "Ring Group"
		if self.rgId != None:
			ret += " (#%d)" % self.rgId
		ret += "\n Rings / Fused To:"
		for i in range(len(self)):
			ret += "\n  * %d %s :\t%s" % (i, repr(self[i]), self.fusedTo[i])
		return ret

class RingSystem(object):
	"""
	Contains all of the rings in the molecule.
	"""

	def __init__(self):
		"""Ring system CTOR."""
		# Smallest Set of Smallest Rings.
		# This may not ultimately be the ring set we draw with.
		#self.sssr = None
		#self.sssrAlgo = None
		# Ring groups -- rings that are touching each other.
		#self.groups = None
		pass

def partition_rings(ringList):
	"""
	Segment a list of rings in the molecule into ring groups based on
	their connectivity. Rings that are not clustered should not be put
	into the same group. 

	Input: A list of Ring objects.
	Output: A list of RingGroup objects.

	This is not based on literature, but aids in analysis and
	construction.
	"""

	def in_same_group(ring1, ring2, ringList):
		"""
		Determine if ring1 and ring2 are in the same group by
		attempting to traverse from one to the other. This is
		just BFS without the queue.
		"""
		# Easy case -- see if the rings are directly connected.
		if ring1.isSpiroTo(ring2) or ring1.isFusedTo(ring2) or \
				ring1.isBridgedTo(ring2):
					return True

		# Iteratively collect rings that are fused/spiro/bridged from
		# ring1. If we manage to collect ring2, then the two are in the
		# same connectivity group. This is similar to [Helson]'s 
		# "central" test algorithm (Procedure D of Ring Analysis, 
		# p. 334)
		rings = list(ringList[:])
		rings.pop(rings.index(ring1))

		collected = [ring1]
		justAdded = [ring1]
		while len(justAdded) != 0:
			remRings = [x for x in rings if x not in collected]
			added = []
			for r in justAdded:
				for s in remRings:
					if r.isSpiroTo(s) or r.isFusedTo(s) or r.isBridgedTo(s):
						collected.append(s)
						added.append(s)
			justAdded = added

		# Result.
		return ring2 in collected

	# Simple case -- only one ring in the entire molecule.
	if len(ringList) == 1:
		return [RingGroup(ringList)]

	# More than one ring in molecule -- must group memberships.
	groups = []
	ungrouped = list(ringList[:])
	justGrouped = []
	rId = 0

	while len(ungrouped) > 0:
		groupRoot = ungrouped.pop(0)
		justGrouped = []
		for ring in ungrouped:
			if in_same_group(groupRoot, ring, ringList):
				justGrouped.append(ring)

		# This is the next ring group.
		group = [groupRoot]
		group.extend(justGrouped)
		groups.append(RingGroup(group, rId))
		rId += 1

		ungrouped = [x for x in ungrouped if x not in justGrouped]

	return groups


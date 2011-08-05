"""
Ring-related data structures.
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

# FIXME: Remove most of Ring's methods. They belong in analysis only!
class Ring(tuple):
	"""
	Ring Data Structure.
	Contains the atom cycle in the ring, bonds, etc.
	Also assists in processing tasks.
	"""

	def __new__(cls, path):
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

	def __init__(self, path):
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
				bonds.append(set([a, a2]))
			return tuple(bonds)

		# XXX/FIXME/TODO: DEPRECATED! This is handled by the tuple superclass.
		self.ringPath = self[:]

		# Bonds in the ring. 
		# A tuple of sets (each set is an edge), which makes testing
		# for bridged rings easy. 
		self.bonds = build_bonds(self.ringPath)

		# The type of ring (elucidated during ring peeling).
		self.type = RING_TYPES.NONE

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

	def isBridgedTo(self, otherRing):
		"""
		Returns True if the ring bridges the other ring.
		"""
		# TODO: Bridged ring support
		#print "TODO: Bridged ring support"
		return False

	def __repr__(self):
		"""Returns representation of object."""
		return str(self)

	def __str__(self):
		"""Returns string representation of the object."""
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

class RingGroup(tuple):
	"""
	Contains all rings that are connected together.
	Eg. "PhPh-C-Ph-C-PhPhPh" contains 3 ring groups.

	As this is a tuple subclass, Pythonic usage should work well:

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
		for ring in self:
			ret += "\n  * " + str(ring)
		return ret

class RingSystem(object):
	"""
	Contains all of the rings in the molecule.
	"""
	pass


"""
Ring analysis.
This code is adapted from [Helson].
"""

from util.enum import enum

RING_TYPES = enum(
	'CORE', 		# Last remaining ring (tractable case)
	'TOUGH_CORE',	# Remaining ring(s) (intractable case)
	'FUSED', 'BRIDGED', 'SPIRO', # Types of peeled rings 
	'IRREGULAR',	# Irregular ring
	'NONE'			# No type assigned
)

class Ring(object):
	"""
	Represents a ring.
	"""

	def __init__(self, path):
		"""Constructor."""

		def build_path(path):
			"""
			Duplicate the list, then remove the end atom if it is also
			the start atom.
			"""
			newPath = path[:]
			if newPath[0] == newPath[-1]:
				newPath.pop()
			return newPath

		def build_bonds(path):
			"""Build a list of sets that represent the bonds."""
			l = len(path)
			bonds = []
			for i in range(l):
				a = path[i]
				a2 = path[(i+1)%l]
				bonds.append(set([a, a2]))
			return bonds

		# The path of atoms in the ring. 
		#	* Each atom only occurs once, so start and end are bonded.
		#	* Technically, rings are undirected cycles.
		self.ringPath = build_path(path)

		# Bonds in the ring. A list of sets; each set is an edge.
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
		rings = ringList[:]
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

	def __str__(self):
		"""Returns string representation of the object."""
		# Type => String map. Don't handle 'none'
		rtype = ''
		if self.type == RING_TYPES.CORE:
			rtype = 'core'
		elif self.type == RING_TYPES.TOUGH_CORE:
			rtype = 'tough_core'
		elif self.type == RING_TYPES.FUSED:
			rtype = 'fused'
		elif self.type == RING_TYPES.BRIDGED:
			rtype = 'bridged'
		elif self.type == RING_TYPES.SPIRO:
			rtype = 'spiro'
		elif self.type == RING_TYPES.IRREGULAR:
			rtype = 'irregular'

		if rtype:
			return "<Ring/%s(%d): %s>"  % (rtype, len(self.ringPath), 
					self.ringPath)

		return "<Ring(%d): %s>" % (len(self.ringPath), self.ringPath)

# ===========================

def ring_analysis(rings):
	"""
	Perform ring analysis using the ring peeling technique.
	This is a two-stage algorithm adapted from [Helson].
	"""
	remainingRings = []
	for r in rings:
		remainingRings.append(Ring(r))

	peelOrder = assign_ring_types(remainingRings)

	# TODO: Repeat again with license for bridged systems. 

	return peelOrder

def assign_ring_types(remainingRings):
	"""
	TODO: DOC
	"""

	def select_fused_spiro_ring(remainingRings):
		"""
		Select the next remaining fused (or spiro) ring to peel.
		The ring chosen has the smallest connectivity with the other
		remaining rings.
		"""
		bestRingPos = None
		minCount = 10000000 # Nothing should have this many shared edges

		for i in range(len(remainingRings)):
			ring = remainingRings[i]

			# Ring cannot be a central ring. 
			if ring.isCentralRing(remainingRings):
				continue

			# Ring cannot be a bridge to other rings. 
			# FIXME: Can speed this up by taking out both rings at once
			isBridge = False
			for r in remainingRings:
				if r is ring:
					continue
				if ring.isBridgedTo(r):
					isBridge = True
					break

			if isBridge:
				continue

			# Count the number of bonds shared with other unassigned 
			# rings.
			# FIXME: Wickedly inefficient. 
			count = 0
			for bond in ring.bonds:		
				for r in remainingRings:
					if ring is r:
						continue
					if bond in r.bonds:
						count += 1

			# XXX/TODO/FIXME:
			# I am unclear about Helson's terminology, but it makes 
			# sense just to skip the rings.
			if count > 3:
				# ring.type = RING_TYPES.IRREGULAR
				continue

			if minCount > count:
				minCount = count
				bestRingPos = i

		# Best ring to peel.
		# Could be 'None'
		return bestRingPos

	def select_bridged_ring(remainingRings):
		"""
		Select the next remaining bridged ring to peel.
		The ring chosen has the smallest connectivity with the other
		remaining rings.
		"""
		bestRingPos = None
		minCount = 10000000 # Nothing should have this many shared edges

		for i in range(len(remainingRings)):
			ring = remainingRings[i]

			# Ring cannot be a central ring
			if ring.isCentralRing(remainingRings):
				continue

			# Count the number of bonds shared with other unassigned
			# rings.
			# FIXME: Wickedly inefficient. 
			count = 0
			for bond in ring.bonds:		
				for r in remainingRings:
					if ring is r:
						continue
					if bond in r.bonds:
						count += 1

			if minCount > count:
				minCount = count
				bestRing = ring

		# Best ring to peel.
		# Could be 'None'
		return bestRingPos

	"""
	Analyze and Peel Rings from the Ring System(s).
	"""

	remainingRings = remainingRings[:]
	peelOrder = []

	while True:
		# If only one more ring exists, we're done. 
		# FIXME: Assume license to redesign irregular rings. 
		if len(remainingRings) == 1:
			ring = remainingRings[0]
			remainingRings[0].type = RING_TYPES.CORE
			peelOrder.append(ring)

		# Select the next best fused/spiro ring, if exists. 
		ringPos = select_fused_spiro_ring(remainingRings)
		if ringPos != None:
			ring = remainingRings.pop(ringPos)

			# Determine if fused or spiro to the remaining rings.
			rtype = None
			for r in remainingRings:
				if ring.isFusedTo(r):
					rtype = RING_TYPES.FUSED
					break
				if ring.isSpiroTo(r):
					rtype = RING_TYPES.SPIRO

			if not rtype:
				# XXX: Error!
				print "Could not determine ring type."
				rtype = RING_TYPES.FUSED

			ring.type = rtype
			peelOrder.append(ring)
			continue

		# Select the next best bridged ring.
		ringPos = select_bridged_ring(remainingRings)
		if ringPos != None:

			ring = remainingRings.pop(ringPos)
			ring.type = RING_TYPES.BRIDGED
			peelOrder.append(ring)
			continue

		# Incomplete assignment... we have to terminate.
		print "Incomplete assignment. Break."
		break

	return peelOrder


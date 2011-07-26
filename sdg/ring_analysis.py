"""
Ring analysis.
This code is adapted from [Helson].
"""

from ring import *

def ring_analysis(rings):
	"""
	Perform ring analysis using the ring peeling technique.
	This is a two-stage algorithm adapted from [Helson].
	"""
	remainingRings = []
	for r in rings:
		remainingRings.append(Ring(r))
	
	# TODO: Should this be external? We need to return RingGroups.
	groups = _segment_into_ring_groups(remainingRings)
	print "Ring Groups: %s" % groups

	peelOrder = _assign_ring_types(remainingRings)

	# TODO: Repeat again with license for bridged systems. 

	return peelOrder

def _segment_into_ring_groups(ringList):
	"""
	Segment the list of rings from the perception stage into ring 
	groups based on the connectivity.

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
		rings = ringList[:]
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
	ungrouped = ringList[:]
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

def _assign_ring_types(remainingRings):
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


"""
Ring analysis.
This code is adapted from [Helson].
"""

from ring import Ring, RingGroup, RING_TYPES

def ring_analysis(ringGroups, mol):
	"""
	Perform ring analysis using the ring peeling technique.
	This is a two-stage algorithm adapted from [Helson].
	"""
	
	for ringGroup in ringGroups:
		peelOrder = _assign_ring_types(ringGroup)

		# TODO: Repeat again with license for bridged systems. 
		ringGroup.peelOrder = peelOrder

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

	rings = []
	for ring in remainingRings:
		rings.append(ring)

	# Can't work with RingGroup directly, convert to list.
	remainingRings = list(remainingRings[:])

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
		# TODO/FIXME: print "Incomplete assignment. Break."
		break

	return peelOrder


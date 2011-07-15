# SDG Preassembly Analysis
# From [Helson 1999]

def ring_perception(graph):
	# Ring perception greatly aids in the assembly phase. 
	# TODO: Need ring perception algorithm
	# Sources are [Balducci 1994] (recommended), [Figueras 1996]
	# Review [Downs 1989] (Found a copy!)
	pass

def build_datastructs(graph):
	# Helson p 326: Recommends building data structures to supplement
	# mere connection tables... 
	pass

def atom_prioritize(graph):
	# TODO: Won't work until ring perception works.
	# This prioritizes the drawing of congested atoms first
	pass


CHAIN_FIXED_ANGLE = 120

def chain_perception(graph):
	"""
	Chain perception information aids in the assembly phase.
	But since chains are not prefabricated units (PFUs), this 
	information may be disregarded during assembly (eg. instances of
	congestion, etc).

		Core Chain Atoms:
			* Acylic (Not partially selected???)
			* Cannot be sp-hybridized
			* Has >= 2 neighbors
			* Has >= 1 acyclic neighbor
			* Has >= 1 acylic beta atom
	"""
	# FIXME: For now, there is no ring support/recognition.
	global CHAIN_FIXED_ANGLE

	# Core chain flags. 
	coreChain = [False for i in range(graph.numAtoms())]
	types = ['EqAS' for i in range(graph.numAtoms())]

	"""Determine which atoms belong to a core chain."""

	for atomNum in range(graph.numAtoms()):
		# FIXME: Assume acyclic
		# FIXME: What does Helson mean, "not partially selected"? p340
		# FIXME: Assume neighbors are acyclic 
		neighbors = graph.getNeighbors(atomNum)

		# Atom cannot be sp-hybridized
		# FIXME: hybridization error handling. 
		if graph.getHybridization(atomNum) in ['sp', 'error']:
			continue
	
		# FIXME: Assume neighbors are acyclic 
		if len(neighbors) < 2:
			continue

		# Determine if it has an acyclic beta atom
		# FIXME: At least one beta atom must be acylic
		hasBeta = False
		queue = neighbors[:]
		while len(queue) > 0:
			at = queue.pop(0)
			for n in graph.getNeighbors(at):
				if n != atomNum:
					hasBeta = True
					break

		if not hasBeta:
			continue

		# Assign as a core chain.
		coreChain[atomNum] = True
		types[atomNum] = 'FxAS'


	"""Step2: reeexamine every FxAS atom to determine if still core."""

	coreChainCopy = coreChain[:]

	for i in range(graph.numAtoms()):
		if types[i] == 'EqAS':
			continue

		neighbors = graph.getNeighbors(i)

		# All Substituents are Heteroatoms
		# XXX: I made an assumption/modification here.
		numCore = 0
		numHetero = 0
		alphaNonCore = len(neighbors)
		for n in neighbors:
			# FIXME/XXX: Assuming core chain not counted!
			if coreChainCopy[n]:
				alphaNonCore -= 1
				numCore += 1
			elif graph.atomTypes[n].upper() not in ['C', 'H']:
				numHetero += 1

		if numHetero > 0 and numHetero >= (alphaNonCore): #- numCore):
			coreChain[i] = False
			types[i] = 'EqAS'
			continue

		# Conjestedness heurisic: num beta atoms > 6
		if len(graph.getBetaAtoms(i)) > 6:
			coreChain[i] = False
			types[i] = 'EqAS'
			continue

		# Degree of substituent heteroatoms.
		# If atom has >= 3 primary heteroatoms, reclassify
		# TODO: Consider building and caching a heteroatom list
		cnt = 0
		for n in neighbors:
			if graph.atomTypes[n].upper() not in ['C', 'H'] and \
					graph.getDegree(n) == 1:
						cnt += 1
		if cnt >= 3:
			coreChain[i] = False
			types[i] = 'EqAS'
			continue

		# Reclassify if has >= 3 substituents with pi systems (double,
		# triple bonds).
		cnt = 0
		for n in neighbors:
			if graph.getDegree(n) in ['sp', 'sp2', 'error']:
				cnt += 1
		if cnt >= 3:
			coreChain[i] = False
			types[i] = 'EqAS'
			continue


	print "================"
	print coreChain
	print types


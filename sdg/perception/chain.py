
class Chain(object):
	"""
	Chain Structure
	Contains metadata corresponding to a single chain. 
	"""
	def __init__(self):
		self.path = [] # Atom numbers in the chain path. 
		self.caps = [] # There are two end caps. 
		self.zigzag = [] # Zigzag pattern along the chain from {L,R}. 

def identify_core_chain(graph):
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
	# FIXME: I will have to develop my own heuristics here. Several
	# needs are listed below. 

	# Core chain flags. 
	coreChainFlags = [False for i in range(graph.numAtoms())]
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
		coreChainFlags[atomNum] = True
		types[atomNum] = 'FxAS'


	"""Step2: reeexamine every FxAS atom to determine if still core."""

	coreChainCopy = coreChainFlags[:]

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
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

		# TODO/FIXME: Need a congestedness heuristic for sitations
		# like this one: CCCCC(CCCC)(CCCC)CCCC. 
		
		# Conjestedness heurisic: num beta atoms > 6
		if len(graph.getBetaAtoms(i)) > 6:
			coreChainFlags[i] = False
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
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

		# Reclassify if has >= 3 substituents with pi systems (double,
		# triple bonds).
		cnt = 0
		for n in neighbors:
			if graph.getDegree(n) in ['sp', 'sp2', 'error']:
				cnt += 1
		if cnt >= 3:
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

	# XXX/FIXME: Assumption: coreChainFlags and types are identical. 
	# I'll return the former for now. 
	#print types
	return coreChainFlags

def identify_chains(graph):
	"""
	We need to take core chain and identify the largest contiguous runs.
	"""
	coreChain = identify_core_chain(graph)

	print "Core Chain:"
	print coreChain

	# Single source LONGEST path.
	for i in range(len(coreChain)):
		if not coreChain[i]:
			continue

		#path = []
		#while
		#n = graph.getNeighbors(i)
		

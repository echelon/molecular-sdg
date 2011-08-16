"""
This module deals with Chain Perception.
Chains are identified and returned. 
"""

# FIXME: This code is _really_ messy.

from algo.path import ShortestPaths
from chain import Chain
	
def identify_chains(mol, rings=None):
	"""
	We need to take core chain and identify the largest contiguous runs.
	"""
	
	def build_connect_mat(mol):
		"""
		Build a connection matrix we can use in Floyd's algorithm.
		'Not connected' is represented by infinity.
		"""
		inf = float('Infinity')
		newMat = [[inf for x in range(mol.size)] for y in range(mol.size)]
		for i in range(mol.size):
			for j in range(mol.size):
				val = mol.connectMat[i][j]
				val = inf if not val else 1
				newMat[i][j] = val
		return newMat

	def get_core_vertices(coreFlags):
		"""Get a list of core vertices from core chain flags."""
		core = []
		for i in range(len(coreFlags)):
			if coreFlags[i]:
				core.append(i)
		return core

	def get_noncore_vertices(coreFlags):
		"""Get a list of non-core vertices from core chain flags."""
		core = []
		for i in range(len(coreFlags)):
			if not coreFlags[i]:
				core.append(i)
		return core

	def disconnect_vertices(matrix, removeVerts = []):
		"""Disconnect the vertices specified from the adj matrix."""
		inf = float('Infinity')
		length = len(matrix)
		for vert in removeVerts:
			for i in range(length):
				matrix[i][vert] = inf
				matrix[vert][i] = inf

	def find_maxweight_verts(shortestPaths):
		"""Find the maximum weight in a ShortestPaths object."""
		inf = float('Infinity')
		DNE = [inf, -inf, None]
		length = shortestPaths.size()
		maxWeight = -1
		maxPath = (-1, -1)
		for i in range(length):
			for j in range(length):
				weight = shortestPaths.getWeight(i, j)
				if weight in DNE:
					continue
				if weight > maxWeight:
					maxWeight = weight
					maxPath = (i, j)

		return maxPath

	# Ring atoms cannot be considered.
	# TODO: Currently, only a lists or tuples is supported as ring input.
	ringAtoms = []
	if rings:
		ratoms = []
		for r in rings:
			ratoms += r
		ringAtoms = list(set(ratoms))

	# Find, identify the "core chain" atoms
	coreFlags = _identify_core_chain(mol, ringAtoms)

	# Copy the molecule's connectivity matrix, then specify that all
	# non-"core chain" atoms are to be removed in the first pass.
	connectMat = build_connect_mat(mol)
	removeAtoms = get_noncore_vertices(coreFlags)

	# Get all of the chains, starting with the longest. 
	chains = []
	while True:
		# Disconnect the specified atoms from the connection matrix
		# This ensures no two chains consist of the same atoms. 
		disconnect_vertices(connectMat, removeAtoms)
	
		# (Re)calculate the shortest paths, and extract the longest chain
		s = ShortestPaths(connectMat)
		verts = find_maxweight_verts(s)
		chain = s.findPath(verts[0], verts[1], include=True)

		# FIXME: I assume this is the only condition that breaks the loop...
		if type(chain) != list:
			break

		# Save the chain, then set up the removal of all atoms that
		# were in it. 
		chains.append(chain)
		removeAtoms = chain[:]
		

	return chains

def _identify_core_chain(mol, ringAtoms = []):
	"""
	Chain perception information aids in the assembly phase.
	But since chains are not prefabricated units (PFUs), this 
	information may be disregarded during assembly (eg. instances of
	congestion, etc).

		Core Chain Atoms:
			* Acylic
			* Cannot be sp-hybridized
			* Has >= 2 neighbors
			* Has >= 1 acyclic neighbor
			* Has >= 1 acylic beta atom
			* (Not partially selected ??)
	"""
	# FIXME: I will have to develop my own heuristics here. Several
	# needs are listed below. 

	# Core chain flags. 
	coreChainFlags = [False for i in range(mol.size)]
	types = ['EqAS' for i in range(mol.size)]

	"""Determine which atoms belong to a core chain."""

	for atomNum in range(mol.size):
		# TODO: What does Helson mean, "not partially selected"? p340
		# Atom cannot be a ring atom. 
		if atomNum in ringAtoms:
			continue

		# Atom cannot be sp-hybridized
		# FIXME: hybridization error handling. 
		if mol.hybridizations[atomNum] in ['sp', 'error']:
			continue
	
		neighbors = mol.alphaAtoms[atomNum]

		# Must have at least two neighbors
		if len(neighbors) < 2:
			continue

		# Must have at least one acyclic neighbor
		hasNeighborAcyclic = False
		for n in neighbors:
			if n not in ringAtoms:
				hasNeighborAcyclic = True
				break

		if not hasNeighborAcyclic:
			continue

		# Must have at least one acyclic beta atom
		hasBetaAcyclic = False
		for b in mol.betaAtoms[atomNum]:
			if n not in ringAtoms:
				hasBetaAcyclic = True
				break

		if not hasBetaAcyclic:
			continue

		# Assign as a core chain.
		coreChainFlags[atomNum] = True
		types[atomNum] = 'FxAS'


	"""Step2: reeexamine every FxAS atom to determine if still core."""

	coreChainCopy = coreChainFlags[:]

	for i in range(mol.size):
		if types[i] == 'EqAS':
			continue

		neighbors = mol.alphaAtoms[i]

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
			elif mol.types[n].upper() not in ['C', 'H']:
				numHetero += 1

		if numHetero > 0 and numHetero >= (alphaNonCore): #- numCore):
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

		# TODO/FIXME: Need a congestedness heuristic for sitations
		# like this one: CCCCC(CCCC)(CCCC)CCCC. 
		
		# Conjestedness heurisic: num beta atoms > 6
		if len(mol.betaAtoms[i]) > 6:
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

		# Degree of substituent heteroatoms.
		# If atom has >= 3 primary heteroatoms, reclassify
		# TODO: Consider building and caching a heteroatom list
		cnt = 0
		for n in neighbors:
			if mol.types[n].upper() not in ['C', 'H'] and \
					mol.degrees[n] == 1:
						cnt += 1
		if cnt >= 3:
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

		# Reclassify if has >= 3 substituents with pi systems (double,
		# triple bonds).
		cnt = 0
		#for n in neighbors:
		#	if mol.hybridizations[n] in ['sp', 'sp2', 'error']:
		#		cnt += 1
		if cnt >= 3:
			coreChainFlags[i] = False
			types[i] = 'EqAS'
			continue

	# XXX/FIXME: Assumption: coreChainFlags and types are identical. 
	# I'll return the former for now. 
	#print types
	return coreChainFlags


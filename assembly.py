"""
Assembly Phase
TODO: Doc
"""

from math import * 

# XXX: Algorithm 6
def substituent_angular_spacing(mol, seedAtom, isHeadAtom=False):
	"""
	Yeilds the optimum angular spacing between the remaining unplaced
	substituents. From [Helson]. 
	"""

	angularDemand = 0 # Area taken up by PFUs around seed atom
	numSubstituent = 0 # Number of non-PFU substituents 

	# TODO TODO TODO: Must be 'global' throughout assembly phase
	# Bonds that have been placed 
	# Bond is 'unplaced' if either of its atoms is unplaced. 
	placedBonds = [] 

	p = [] # Set of adjacent bonds (in PFUs) that have been encountered

	# Work with adjacent, unplaced bonds.
	bonds = [frozenset((n, seedAtom)) for n in mol.alphaAtoms[seedAtom]]
	bonds = filter(lambda b: b not in placedBonds, bonds)
	
	for bond in bonds:

		print bond

		if bond in p:
			continue

		PFU = [] # XXX TEMP
		if bond in PFU:
			g = iter(bond - set([seedAtom])).next()
			# TODO:
			# Increase AngleDemand by demand of PFU, which is the
			# complement of the CFS at atom g
			
			# TODO: Add all bonds in PFU to set P

			pass
			
		else:
			numSubstituent += 1

		
	if isHeadAtom:
		numSubstituent -= 1

	# TODO
	# Correct for FxAS. If the seed atom is a core chain atom,
	# decrement NumSub and increase angular demabd by chain angle.


	# TODO
	# angularSpacing = (seedAtom.cfs - angularDemand) / (numSub + 1)

# XXX: Algorithm 7
def substituent_sequence(mol, seedAtom):

	remaining = []
	complete = []

	# Return sequence
	sequence = []

	# TODO: Initialize to the set of atoms that have and have not been placed.
	remaining = list(mol.alphaAtoms[seedAtom])

	while len(remaining) > 0:

		s = None

		# TODO
		chain = False # XXX TEMP
		if seedAtom is chain:
			pass # TODO

		if len(complete) == 0:
			pass # TODO

		else:
			#s = substituent immediately CCW of last chosen atom
			pass # TODO

		"""
		ensure that s is not from the wrong side of a PFU --
			let b be the bond from the seed atom to s.
			if b is in a PFU:
				examine the PFU's bonds that are adj to the seed atom,
				and locate the bond q that is equal to the seed atom's local
				CFS.hi (there must be one).
				set s to q's other atom. 
				add the PFU's atoms to complete and subtract the PFU's 
				atoms from  remaining
		"""

		complete.append(s)
		remaining.remove(s)
		sequence.append(s)

	return sequence

# XXX: Algorithm 8
def substituent_placement_for_ring_atom_seed():
	pass

# XXX: Algorithm 9
def substituent_placement_for_core_chain_atom_seed():
	pass

# XXX: Algorithm 10
def substituent_placement_for_uncategorized_seed_atom():
	pass

# XXX: Algorithm 11
def place_atom():
	pass
	

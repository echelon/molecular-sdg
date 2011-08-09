#!/usr/bin/env python

"""
This represents the current work on the Molecular Graph and Structure
Diagram Generation (SDG) work that I am doing. There are several files
in this repository that are not relevant anymore, but contain an earlier
approach at drawing that I wish to retain for the time being.

At present, the code I am working on is in 'matrix.py' and 'smiles.py'.
"""

import sys

from molecule import Molecule
from ring import partition_rings

from examples import get_example
from smiles import Smiles
from smiles import smiles_to_molecule
from algo.path import *

# Perception
from perception.rings import *
from perception.chain import *

# Analysis Phase
from sdg.ring_analysis import *
from sdg.ring_construction import *

# XXX: Temp
from gui import Window

def main():
	"""Main function"""

	smiles = None
	name = None
	if len(sys.argv) < 2:
		ex = get_example()
		name = ex[1]
		smiles = ex[2]
		print ">>> Need to supply SMILES text as argument."
		print ">>> Using %s \"%s\" as an example.\n" % ex[1:] 
	else:
		ex = get_example(sys.argv[1])
		if ex:
			name = ex[1]
			smiles = ex[2]
			print ">>> Using %s per argument.\n" % name
		else:
			smiles = sys.argv[1]

	mol = smiles_to_molecule(smiles)
	mol.informalName = name

	#mol.print_matrix()

	# Perception algorithms. 
	rings = identify_rings(mol)
	chains = identify_chains(mol, rings)

	# Ring analysis
	ringGroups = partition_rings(rings)

	#print ringGroups

	ring_analysis(ringGroups, mol)

	for rg in ringGroups:
		print rg.peelOrder

	for group in ringGroups:
		print group
		construct_group(group)

	win = Window('title')
	win.setText(smiles)
	win.ringGroups = ringGroups # XXX: Pass ring groups off to gui
	win.run()


if __name__ == '__main__':
	main()


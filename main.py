#!/usr/bin/env python

"""
This represents the current work on the Molecular Graph and Structure
Diagram Generation (SDG) work that I am doing. There are several files
in this repository that are not relevant anymore, but contain an earlier
approach at drawing that I wish to retain for the time being.

At present, the code I am working on is in 'matrix.py' and 'smiles.py'.
"""

# Std lib
import sys

# Project
from examples import get_example
from molecule import Molecule
from smiles import Smiles
from smiles import smiles_to_molecule
from util.matrix import print_matrix
from algo.path import *
from perception.rings import *
from perception.chain import *
from sdg.ring_analysis import *

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

	mol.print_matrix()

	# Perception algorithms. 
	rings = identify_rings(mol)
	chains = identify_chains(mol, rings)

	print "\nRing Perception:"
	print rings
	print "\nChain Perception:"
	print chains

	# Ring analysis
	peelOrder = ring_analysis(rings, mol)	

	print peelOrder

	sys.exit()

	# TODO: Remove 'Molecule' 
	#print "\n\nRing Peel Order:"
	#for x in peelOrder:
	#	print x

	win = Window('title')
	win.setText(smiles)
	win.run()

	#print "\n\n(Repeated) Ring Perception (%d):\n%s\n" % (len(rings), rings)
	#print "Smiles: %s" % smiles

if __name__ == '__main__':
	main()


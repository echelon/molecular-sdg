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
import random

# Project
from matrix import MolMatrix
from smiles import Smiles
from smiles import smiles_to_matrix
from sdg.analysis import *

def get_example():
	"""Return an example molecule (name, smiles) tuple."""

	EXAMPLES = {
		'hexane': 'CCCCCC',
		'hexene': 'C=CCCCC',
		'hexyne': 'C#CCCCC',
		'isohexane': 'CC(C)CCC',
		'adenine': 'n1c(c2c(nc1)ncn2)N',
		'adenine{2}': 'c1[nH]c2c(ncnc2n1)N', # XXX: Error! Hydrogen included.
		'acetic acid': 'CC(=O)O',
		'p-bromocholobenzene': 'C1=CC(=CC=C1Cl)Br',
		'vanillin': 'O=Cc1ccc(O)c(OC)c1',
		'naphthalene': 'c1cccc2c1cccc2',
	}

	key = random.choice(EXAMPLES.keys())
	return (key, EXAMPLES[key])

def main():
	"""Main function"""

	smiles = None
	if len(sys.argv) < 2:
		ex = get_example()
		print ">>> Need to supply SMILES text as argument."
		print ">>> Using %s \"%s\" as an example.\n" % ex 
		smiles = ex[1]
	else:
		smiles = sys.argv[1]

	mol1 = smiles_to_matrix(smiles)
	mol1.print_matrix()
	print "\n"

	#mol2 = mol1.canonicalize()
	#mol2.print_matrix()
	#print "\n"

	# XXX: Testing...
	print "Chain Perception: \n"
	chain_perception(mol1)
	print smiles

if __name__ == '__main__':
	main()


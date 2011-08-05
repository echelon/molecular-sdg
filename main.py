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
from matrix import MolMatrix, ConstMolMatrix
from smiles import Smiles
from smiles import smiles_to_matrix
from util.matrix import print_matrix
from algo.path import *
from perception.rings import *
from perception.chain import *
from sdg.ring_analysis import *

# XXX: Temp
from gui import Window

def get_example(key=None):
	"""
	Return an example molecule (name, smiles) tuple.
	Can provide a lookup key, or get a random molecule.
	"""

	EXAMPLES = {
		'hexane': 'CCCCCC',
		'hexene': 'C=CCCCC',
		'hexyne': 'C#CCCCC',
		'isohexane': 'CC(C)CCC',
		'adenine': 'n1c(c2c(nc1)ncn2)N',
		'adenine2': 'c1[nH]c2c(ncnc2n1)N', # XXX: Hydrogen included!
		'acetic acid': 'CC(=O)O',
		'p-bromocholobenzene': 'C1=CC(=CC=C1Cl)Br',
		'vanillin': 'O=Cc1ccc(O)c(OC)c1',
		'toulene': 'Cc1ccccc1',

		# Polycyclic Aromatic Hydrocarbons
		'naphthalene': 'c1cccc2c1cccc2',
		'acenaphthene': 'c2cc1cccc3c1c(c2)CC3',
		'benzo[k]fluoranthene': 'c1ccc2cc-3c(cc2c1)-c4cccc5c4c3ccc5',
		'dibenz(a,h)anthracene': 'c1ccc2c(c1)ccc3c2cc4ccc5ccccc5c4c3',
		'coronene': 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67', # 'superbenzene'
		'dicoronylene': 'c1cc2ccc3cc4c9cc%10ccc%11ccc%12ccc%13ccc%14cc(c5cc6'
				+ 'ccc7ccc1c8c2c3c(c45)c6c78)c9c%15c%10c%11c%12c%13c%14%15',

		# Basic Stereochemistry
		'coenzyme-a': 'c1c2c(cc(c1F)N3CCNCC3)n(cc(c2=O)C(=O)O)C4CC4',
		'coenzyme-a2': 'O=C(NCCS)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC'
				+ '[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3OP(=O)(O)O',

		# Bridged systems
		'pericine': 'c23c1ccccc1nc2\C(=C)[C@H]4C(=C/C)\CN(CC3)CC4',
		'trogers-base': 'c1(ccc3c(c1)CN4c2ccc(cc2CN3C4)C)C',
		'lurasidone': 'O=C1N(C(=O)[C@H]3[C@@H]1[C@@H]2CC[C@H]3C2)C[C@@H'
						+ ']4CCCC[C@H]4CN7CCN(c6nsc5ccccc56)CC7',

	}

	def get_key(k, li):
		# Random or exact key
		if key == None:
			return random.choice(li.keys())
		if key in li:
			return key

		# Key is a numeric offset into list 
		# XXX: Bad idea! These aren't stored in order.
		try:
			if str(int(k)) == k:
				return li.keys()[int(k)]
		except:
			pass

		return None

	key = get_key(key, EXAMPLES)

	if not key:
		return None

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
		ex = get_example(sys.argv[1])
		if ex:
			print ">>> Using %s per argument.\n" % ex[0] 
			smiles = ex[1]
		else:
			smiles = sys.argv[1]

	mol = smiles_to_matrix(smiles)

	mol = ConstMolMatrix(mol)
	mol.print_matrix()

	# Perception algorithms. 
	rings = identify_rings(mol)

	print rings
	for ring in rings:
		print ring

	chains = identify_chains(mol, rings)

	#print "SMILES: %s\n" % smiles
	#print "Ring Perception (%d):\n%s\n" % (len(rings), rings)
	#print "\nChain Perception (%d):\n%s\n" % (len(chains), chains)

	# Ring analysis
	peelOrder = ring_analysis(rings, mol)	

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


"""
List of example molecules to test, as well as a random selection
procedure.
"""

import random

# Example molecules, grouped by class/complexity. 
EXAMPLES = {
	# Basic alkanes, alkenes, alkynes and branched hydrocarbons
	'hydrocarbons': {
		'hexane': 'CCCCCC',
		'hexene': 'C=CCCCC',
		'hexyne': 'C#CCCCC',
		'isohexane': 'CC(C)CCC',
	},

	# Basic organic chemistry
	'organic': {
		'acetic acid': 'CC(=O)O',
		'p-bromocholobenzene': 'C1=CC(=CC=C1Cl)Br',
		'vanillin': 'O=Cc1ccc(O)c(OC)c1',
		'toulene': 'Cc1ccccc1',
	},

	# Basic bio molecules
	'biochem': {
		'adenine': 'n1c(c2c(nc1)ncn2)N',
		'adenine2': 'c1[nH]c2c(ncnc2n1)N', # XXX: Hydrogen included!
	},

	# Polycyclic Aromatic Hydrocarbons
	'pah': {
		'naphthalene': 'c1cccc2c1cccc2',
		'acenaphthene': 'c2cc1cccc3c1c(c2)CC3',
		'benzo[k]fluoranthene': 'c1ccc2cc-3c(cc2c1)-c4cccc5c4c3ccc5',
		'dibenz(a,h)anthracene': 'c1ccc2c(c1)ccc3c2cc4ccc5ccccc5c4c3',
		'coronene': 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67', # 'superbenzene'
		'dicoronylene': 'c1cc2ccc3cc4c9cc%10ccc%11ccc%12ccc%13ccc%14cc(c5cc6'
				+ 'ccc7ccc1c8c2c3c(c45)c6c78)c9c%15c%10c%11c%12c%13c%14%15',
	},

	# Exhibit stereochemistry
	'sterochem' : {
		# bio molecules with stereochemistry
		'coenzyme-a': 'c1c2c(cc(c1F)N3CCNCC3)n(cc(c2=O)C(=O)O)C4CC4',
		'coenzyme-a2': 'O=C(NCCS)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC'
				+ '[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3OP(=O)(O)O',
	},

	# Bridged systems. Most also exhibit sterochemistry. Very complex!
	'bridged': {
		'pericine': 'c23c1ccccc1nc2\C(=C)[C@H]4C(=C/C)\CN(CC3)CC4',
		'trogers-base': 'c1(ccc3c(c1)CN4c2ccc(cc2CN3C4)C)C',
		'lurasidone': 'O=C1N(C(=O)[C@H]3[C@@H]1[C@@H]2CC[C@H]3C2)C[C@@H'
						+ ']4CCCC[C@H]4CN7CCN(c6nsc5ccccc56)CC7',
	},
}


def get_example(key=None):
	"""
	Return an example molecule (name, smiles) tuple.
	Can provide a lookup key, or get a random molecule.
	"""

	k1 = None # Outer key
	k2 = None # Inner key
	val = None # SMILES value

	# If outer key specified, get random inner value. 
	if key in EXAMPLES:
		k1 = key
		k2 = random.choice(EXAMPLES[key].keys())
		val = EXAMPLES[k1][k2]
		return (k1, k2, val)

	innerKeys = {}
	for outk in EXAMPLES:
		for ink in EXAMPLES[outk]:
			innerKeys[ink] = outk

	# If inner key specified, get exact inner value.
	if key in innerKeys:
		k1 = innerKeys[key]
		k2 = key
		val = EXAMPLES[k1][k2]
		return (k1, k2, val)
	
	# Return a random value.
	key = random.choice(innerKeys.keys())
	k1 = innerKeys[key]
	k2 = key
	val = EXAMPLES[k1][k2]
	return (k1, k2, val)


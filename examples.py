"""
List of example molecules to test, as well as a random selection
procedure.
"""

import random

# Example molecules, grouped by class/complexity.
# The 'names' used as keys are not necessarily official nomenclature.
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
		'benzenesulfonic acid': 'c1ccc(cc1)S(=O)(=O)O',
		'o-xylene': 'Cc1ccccc1C',
		'm-xylene': 'Cc1cccc(c1)C',
		'p-xylene': 'Cc1ccc(cc1)C',
		'dinitroaniline': 'O=[N+]([O-])c1cc(ccc1N)[N+]([O-])=O', 
	},

	# Basic bio molecules
	'biochem': {
		'adenine': 'n1c(c2c(nc1)ncn2)N',
		'adenine2': 'c1[nH]c2c(ncnc2n1)N', # XXX: Hydrogen included!
		'guanine': 'c1[nH]c2c(n1)c(=O)[nH]c(n2)N',
		'cytosine': 'c1cnc(=O)[nH]c1N',
		'thymine': 'Cc1c[nH]c(=O)[nH]c1=O',
		'uracil': 'c1c[nH]c(=O)[nH]c1=O',
	},

	# Polycyclic Aromatic Hydrocarbons
	'pah': {
		'naphthalene': 'c1cccc2c1cccc2', # 2
		'acenaphthene': 'c2cc1cccc3c1c(c2)CC3', # 3
		'chrysene': 'c1ccc2c(c1)ccc3c2ccc4c3cccc4', # 4
		'dibenz(a,h)anthracene': 'c1ccc2c(c1)ccc3c2cc4ccc5ccccc5c4c3', # 5
		'benzo[k]fluoranthene': 'c1ccc2cc-3c(cc2c1)-c4cccc5c4c3ccc5', # 5
		'coronene': 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67', # 'superbenzene'
		'dicoronylene': 'c1cc2ccc3cc4c9cc%10ccc%11ccc%12ccc%13ccc%14cc(c5cc6'
				+ 'ccc7ccc1c8c2c3c(c45)c6c78)c9c%15c%10c%11c%12c%13c%14%15',
	},
	
	# Multiple ring groups
	'multiring': {
		'...trichloro-dibenzene': 'ClC(Cl)(Cl)C(c1ccccc1)c2ccccc2',
		'...tetrahydronaphthalene': 'c3ccc2c(CCCC2CCc1ccccc1)c3', # 2, 1 rings
		'...tetrahydroisoquinolinium': 'O=C(OCCCCCOC(=O)CC[N+]2(C(c1c(cc(OC)'
			+ 'c(OC)c1)CC2)Cc3ccc(OC)c(OC)c3)C)CC[N+]5(C)C(c4cc(OC)c(OC)'
			+ 'cc4CC5)Cc6ccc(OC)c(OC)c6', # 2, 2, 1, 1 rings
		'...azo-benzenesulfonic acid': 'c1ccc2c(c1)ccc(c2N=Nc3ccc(cc3S(=O)'
			+ '(=O)O)N=Nc4ccc(cc4)S(=O)(=O)O)O', # 2, 1, 1  rings
	},

	# Spiro compounds
	'spiro': {
		'2-azaspiro[4.5]decan-3-one': 'O=C2NCC1(CCCCC1)C2',
		# Fluorescein is conjested!!
		'fluorescein': 'c1ccc2c(c1)C(=O)OC23c4ccc(cc4Oc5c3ccc(c5)O)O',
		'irbesartan': 'O=C1N(\C(=N/C12CCCC2)CCCC)Cc5ccc(c3ccccc3c4nnnn4)cc5',
	},

	# Exhibit stereochemistry
	'sterochem' : {
		# bio molecules with stereochemistry
		'coenzyme-a': 'c1c2c(cc(c1F)N3CCNCC3)n(cc(c2=O)C(=O)O)C4CC4',
		'coenzyme-a2': 'O=C(NCCS)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC'
				+ '[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3OP(=O)(O)O',
		'atp': 'c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO'
			+ '[P@@](=O)(O)O[P@@](=O)(O)OP(=O)(O)O)O)O)N',
		'luciferin': 'O=C(O)[C@@H]1NC(/SC1)=C2/S\C\3=C\C(=O)\C=C/C/3=N2',
		# Drug with spiro connectivity
		'griseofulvin':'O=C2c3c(O[C@@]21C(/OC)=C\C(=O)C[C@H]1C)c(Cl)c(OC)cc3OC',
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


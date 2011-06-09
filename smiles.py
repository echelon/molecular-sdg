#

def tokenize_smiles(smileStr):
	"""
	Convert SMILES string into a list of tokens. 
	Important because 'Cl', 'Br', etc. constitute one token.
	"""
	READ_AHEAD = {
		'C': 'Cl',
		'B': 'Br'
	}

	queue = []
	pos = 0
	while pos < len(smileStr):
		ch = smileStr[pos]

		# Check for Br, Cl, etc.
		if ch in READ_AHEAD and pos < len(smileStr)-1:
			rd = ch + smileStr[pos+1]
			if rd == READ_AHEAD[ch]:
				ch = rd
				pos += 1

		# TODO: Connectivity beyond '9'

		queue.append(ch)
		pos += 1

	return queue

class Smiles(object):
	"""
	Work with a SMILES string.
	TODO: Canonicalization.
	"""

	def __init__(self, smilesStr):
		"""CTOR."""

		self.string = smilesStr
		self.tokens = tokenize_smiles(smilesStr)
		self.canonicalized = None
		self.canonicalized_tokens = None


	def num_atoms(self):
		"""Calculate the number of non-H atoms"""
		cnt = 0
		for x in self.tokens:
			x = str(x)
			if x.isalpha() and x.upper() != 'h':
				cnt += 1

		return cnt


from matrix import *
import string
import sys # TODO: Temp for sys.exit()

class Smiles(object):
	"""
	SMILES (Simplified Molecular Input Line Entry Specification)

	This class parses SMILES strings into a collection of tokens that
	are easily converted into molecular graphs.
	"""

	def __init__(self, smilesStr):
		"""CTOR."""

		self.string = smilesStr
		self.tokens = self.tokenizeString(smilesStr)

		# TODO: Canonoicalization
		# XXX: Perhaps this isn't something I should even manage here...
		self._cachedCanonical = None # Should be 'Smiles' instance

	def numAtoms(self):
		"""Calculates the number of non-H atoms"""
		# XXX: This hydrogen-excluded policy WILL cause problems in the 
		# future if I don't develop a comprehensive protocol for
		# handling.
		cnt = 0
		for x in self.tokens:
			x = str(x)
			if x.isalpha() and x.upper() != 'h':
				cnt += 1
		return cnt

	def __str__(self):
		"""String Representation."""
		return "<Smiles: %s>" % self.string

	@staticmethod
	def tokenizeString(smileStr):
		"""
		Convert SMILES string into a list of tokens. Dealing with
		tokens is much more convenient than working with a character
		array. Nothing is removed, simply grouped.

		The following items are multiple characters wide, but constitute
		one token here:

			* Organic subset atoms 'Cl' and 'Br'
			* Digits beyond 10 (TODO)
			* Charges (various formats): -, +++, 3-, +2
		
		Hydrogen and Inorganic atoms require special processing here.
		(TODO)
		"""
		# TODO: This won't work for inorganic set (bracketed) 
		# TODO: Stereochemistry
		# TODO: This won't work for ring closures beyond the 9th.
		READ_AHEAD = {
			'C': 'Cl',
			'B': 'Br'
		}
		CHARGE = '+-' + string.digits
		CHARGE_SIGN = '+-'

		print smileStr

		# Process string. 
		tokens = []
		pos = 0
		inBracket = False
		while pos < len(smileStr):
			char = smileStr[pos]

			# Current bracket state
			if char == '[':
				inBracket = True
				tokens.append(char)
				pos += 1
				continue

			if char == ']':
				inBracket = False
				tokens.append(char)
				pos += 1
				continue

			# Group atom charges.
			# Charges only occur in brackets. (TODO: Confirm via spec.)
			if char in CHARGE and inBracket:
				# Only confirm that this is a charge when a '+' or '-' 
				# is encountered. Otherwise, digits could mean anything.
				# XXX: This algorithm allows for invalid formats such 
				# as '++5++2' and so forth, but at present I don't care
				isCharge = False
				chargeStr = ''
				for i in range(pos, len(smileStr)):
					rd = smileStr[i]
					if rd not in CHARGE:
						break
					elif rd in CHARGE_SIGN:
						isCharge = True

					chargeStr += rd

				if isCharge:
					tokens.append(chargeStr)
					pos += len(chargeStr)
					continue

			# Handle two+ character Inorganic set (eg. Co, Sn, etc.)
			# We must be VERY careful not to group a bonded hydrogen.
			# Examples of special cases to worry about are [OH] 
			# (Hydroxide) and [nH] (Aromatic Nitrogen bonded to a  
			# Hydrogen). Note: Inorganic only occur in brackets!
			if char in string.uppercase and inBracket:
				# TODO: Document design rationale
				isInorganic = True
				atomStr = ''
				for i in range(pos, len(smileStr)):
					rd = smileStr[i]
					if rd not in string.letters:
						break
					atomStr += rd

				if len(atomStr) > 1:
					# Already know first letter is uppercase, but
					# all subsequent letters must be lowercase
					for x in atomStr[1:]:
						if x not in string.lowercase:
							isInorganic = False

					if isInorganic:
						tokens.append(atomStr)
						pos += len(atomStr)
						continue

			# Handle Br and Cl (Organic Subset)
			if char in READ_AHEAD and pos < len(smileStr)-1:
				rd = char + smileStr[pos+1]
				if rd == READ_AHEAD[char]:
					tokens.append(rd)
					pos += 2
					continue

			# Catch other cases
			tokens.append(char)
			pos += 1

		print tokens


		#print "Exiting"
		#sys.exit() # XXX: TEMPORARY
		return tokens 

	def toMatrix(self):
		"""Convert a SMILES string into a adjacency matrix representation."""

		# First, we must convert the input string into a proper queue
		# Organic subset: B, C, N, O, P, S, F, Cl, Br, I 
		# Everything else, including H: must be specified in brackets.
		ATOMS = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
		ATOMS_UPPER = map(lambda x: x.upper(), ATOMS)

		# Additional Symbols
		BRANCHING = '()' # Branch Start & End
		BONDS = '=#'	 # Double and Triple Bonds
		CONNECTIVE = '%' # Connectivity beyond '9' -- TODO NOT YET HANDLED.

		queue = self.tokens
		mat = MolMatrix(self.numAtoms())

		# Symbol type tests
		def is_atom(sym): return sym.upper() in ATOMS_UPPER
		def is_bond(sym): return sym in BONDS
		def is_branch_start(sym): return sym == '('
		def is_branch_end(sym): return sym == ')'
		def is_closure(sym): return sym.isdigit()
		def is_bond_order(sym): return sym in '=#'

		# Make note of the connection beween two atoms. 
		def connect(a1, a2, bondOrder=1):
			mat.connectMat[a1][a2] = True
			mat.connectMat[a2][a1] = True
			mat.bondOrderMat[a1][a2] = bondOrder
			mat.bondOrderMat[a2][a1] = bondOrder

		# Label the atom type
		def label(a, name):
			mat.atomTypes[a] = name

		# Keep track of state during parsing
		atomCnt = 0 # Total num atoms
		atomPrev = 0 # Previous atom 
		bondOrder = 1 # Bond order: 1, 2, 3
		branchStack = [] # Branching, eg. C(C)(C)C
		ringClosures = {} # Keep track of cycles, eg. c1ccccc1

		for sym in queue:

			if is_atom(sym):
				atomId = atomCnt
				atomCnt += 1
				label(atomId, sym)

				if atomId != 0:
					# TODO: Conjugated systems bond order = 1.5
					connect(atomPrev, atomId, bondOrder)
					bondOrder = 1

				atomPrev = atomId

			elif is_branch_start(sym):
				branchStack.append(atomPrev)

			elif is_branch_end(sym): 
				atomPrev = branchStack.pop()

			elif is_closure(sym):
				# TODO: Won't handle > 9
				num = int(sym)
				if num not in ringClosures:
					ringClosures[num] = atomPrev
				else:
					# TODO: Conjugated systems bond order = 1.5
					connect(atomPrev, ringClosures[num], bondOrder)
					bondOrder = 1

			elif is_bond_order(sym):
				if sym == '=':
					bondOrder = 2
				else:
					bondOrder = 3

		# End queue processing, return matrix.
		return mat

def smiles_to_matrix(smiles):
	"""Function to return a MolMatrix from a smiles string without an
	intermediate."""
	# XXX: Do I really need to keep this?

	if type(smiles) == str:
		smiles = Smiles(smiles)

	return smiles.toMatrix()

def matrix_to_smiles(matrix):
	# TODO: Need to convert back into a smiles expression.
	pass


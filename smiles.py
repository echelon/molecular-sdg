from matrix import *
import string

# TODO: Reorganize class
# TODO: Fix documentation 
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
		Convert SMILES string into a list of tokens. Working with
		tokens is much more convenient than dealing with a character
		array. Nothing is removed from the string, simply grouped.

		The following items are multiple characters wide, but constitute
		'one token' each:

			* Unbracketed Organic set atoms 'Cl' and 'Br'
			* Bracketed Inorganics, such as 'Co' and 'Sn'
			* Atomic charges (various formats): -, +++, 3-, +2
			* Tetrahedral sterochemistry direction notation '@@'
			* Digits beyond 9
		"""
		READ_AHEAD = {
			'C': 'Cl',
			'B': 'Br'
		}
		CHARGE = '+-' + string.digits
		CHARGE_SIGN = '+-'

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

			# Handle tetrahedral stereochemistry symbol @@
			# These only occur in brackets (TODO: Confirm via spec)
			if char is '@' and inBracket and pos < len(smileStr)-1:
				rd = char + smileStr[pos+1]
				if rd == '@@':
					tokens.append(rd)
					pos += 2
					continue

			# Handle Br and Cl (Organic Subset)
			if char in READ_AHEAD and pos < len(smileStr)-1:
				rd = char + smileStr[pos+1]
				if rd == READ_AHEAD[char]:
					tokens.append(rd)
					pos += 2
					continue

			# Group digits greater than nine.
			# XXX: This collides with some other issues, like charges.
			# Be sure those issues are tokenized first!
			if char in string.digits:
				digitStr = ''
				for i in range(pos, len(smileStr)):
					rd = smileStr[i]
					if rd not in string.digits:
						break
					digitStr += rd

				tokens.append(digitStr)
				pos += len(digitStr)
				continue

			# Catch other cases
			tokens.append(char)
			pos += 1

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
		CHARGE = '+-'    # Cation/anion charge. Only occur in brackets.

		queue = self.tokens
		mat = MolMatrix(self.numAtoms())

		# Symbol type tests
		#def is_atom(sym): return sym.upper() in ATOMS_UPPER
		def is_atom(sym): return type(sym) is str and sym.isalpha() 
		def is_bond(sym): return sym in BONDS
		def is_closure(sym): return sym.isdigit() # TODO: Not a good test
		def is_bond_order(sym): return sym in '=#'
		def is_charge(sym): return len(filter(lambda x: x in CHARGE, sym)) > 0

		# Make note of the connection beween two atoms. 
		def connect(a1, a2, bondOrder=1):
			mat.connectMat[a1][a2] = True
			mat.connectMat[a2][a1] = True
			mat.bondOrderMat[a1][a2] = bondOrder
			mat.bondOrderMat[a2][a1] = bondOrder

		# Label the atom type
		def label(a, name, isotope = 0):
			mat.atomTypes[a] = name
			mat.atomIsotopes[a] = isotope

		# Parse the charge. Examples: +, --, 2+, -3
		# TODO: Make more compact, and more valid to spec.
		def parse_charge(chargeStr):
			charge = 0
			hasDigit = False
			for x in chargeStr:
				if x.isdigit():
					hasDigit = True
					break

			# Type 1: 3+, -2
			if hasDigit:
				charge = int(chargeStr.replace('+', '').replace('-', ''))
				if '-' in chargeStr:
					charge *= -1
				return charge

			# Type 2: +++, --
			for x in chargeStr:
				if x is '-':
					charge -= 1
				elif x is '+':
					charge += 1
			return charge

		# Keep track of state during parsing
		atomCnt = 0 # Total num atoms
		atomPrev = 0 # Previous atom 
		bondOrder = 1 # Bond order: 1, 2, 3
		isotope = 0 # Atom isotope
		branchStack = [] # Branching, eg. C(C)(C)C
		ringClosures = {} # Keep track of cycles, eg. c1ccccc1

		# Bracket state. Brackets denote isotope, charge, inorganics...
		inBrackets = False # Currently in brackets
		inBracketsAtomFound = False # Atom encountered in brackets

		for sym in queue:

			if is_atom(sym):
				atomId = atomCnt
				atomCnt += 1

				label(atomId, sym, isotope)
				isotope = 0

				if inBrackets:
					# Once atom found, isotope can't be set.
					inBracketsAtomFound = True

				if atomId != 0:
					# TODO: Conjugated systems bond order = 1.5
					connect(atomPrev, atomId, bondOrder)
					bondOrder = 1

				atomPrev = atomId

			# Handle bracket state 
			if sym == '[':
				inBrackets = True
				inBracketsAtomFound = False
				continue
			if sym == ']':
				inBrackets = False
				inBracketsAtomFound = False
				continue

			# Handle branching
			if sym == '(':
				branchStack.append(atomPrev)
				continue
			if sym == ')':
				atomPrev = branchStack.pop()
				continue

			if is_closure(sym):
				# TODO: Won't handle > 9
				num = int(sym)
				if num not in ringClosures:
					ringClosures[num] = atomPrev
				else:
					# TODO: Conjugated systems bond order = 1.5
					connect(atomPrev, ringClosures[num], bondOrder)
					bondOrder = 1

			if is_bond_order(sym):
				if sym == '=':
					bondOrder = 2
				else:
					bondOrder = 3
				continue

			if is_charge(sym) and inBrackets:
				# The '-' symbol is also used for a few limited cases
				# of single bond representation. 
				mat.atomCharges[atomId] = parse_charge(sym)
				continue

			# Isotope can only be set in brackets, before the atom.
			if sym.isdigit() and inBrackets and not inBracketsAtomFound:
				isotope = sym

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


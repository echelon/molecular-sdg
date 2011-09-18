"""
TODO: Doc.
"""

from weights import ATOMIC_WEIGHTS
from point import Point

# TODO: Add a __key__ check for assignment
# 		Any invalid keys should not support assignment!
#
# TODO: Add a molecular weight generator. This entails calculating the
#		number of hydrogens. We'll need an 'atoms' dictionary with data.
class Molecule(object):
	"""
	TODO: Docs.
	"""

	def __init__(self, types, bondOrderMat, connectMat=None, charges=None,
				 isotopes=None, ringSystem=None, smiles=None):
		"""
		Molecule constructor.
		Supply input necessary to build the molecule object:

		Mandatory:
			* types	-- Labels of the atoms. C, N, O, etc.
			* bondOrderMat -- Weighted adjacency matrix; weights are
							  the bond orders between atom pairs.
							  Weights : 1, 1.5 (aromatic), 2, and 3. 
		Optional:
			* connectMat -- Boolean adjacency matrix.
			* charges -- Charges on the atoms. Defaults to 0.
			* isotopes -- Atom isotopes. Default to 0, meaning regular.
			* ringSystem -- Rings in the system.
			* smiles -- Smiles text.

		Everything else is calculated from the supplied information.
		"""

		# Number of atoms. Typically non-Hydrogen included.
		self.size = len(connectMat)

		# Connectivity matrix 
		self.connectMat = None

		# Bond orders between atom pairs: 1, 1.5 (aromatic), 2, and 3
		self.bondOrderMat = None

		# Atom labels (tuples)
		self.types = None # C, O, N, Cl, etc.
		self.charges = None
		self.isotopes = None

		# Calculated atom labels
		self.degrees = None
		self.hybridizations = None # sp, sp2, sp3, or 'error'

		# Calculated neighbor and neighbor-of-neighbor tables. 
		self.alphaAtoms = None
		self.betaAtoms = None

		# Rings and chains for the molecule
		self.chains = None
		self.rings = None
		self.ringGroups = None

		# Per-atom flags set in perception, etc. phases.
		self.isInRing = [False for x in range(self.size)]
		self.isInChain = [False for x in range(self.size)] # Not capping subst!

		# Per-atom reference to ring group and chain membership, if exists
		self.chainRef = [None for x in range(self.size)]
		self.ringGroupRef = [None for x in range(self.size)]

		# Reference to smiles text, etc. (optional)
		self.smiles = None
		self.informalName = None

		# Circular Free Sweep (CFS) for each atom. Used in ring
		# construction of analysis phase and later again in assembly.
		self.cfsInitialized = [False for x in range(self.size)]
		self.cfs = [{'hi':0, 'lo':0} for x in range(self.size)]

		# Drawing data --
		self.isPlaced = [False for x in range(self.size)]
		self.pos = [Point() for x in range(self.size)]

		# Calculated number of hydrogens for each atom.
		self.hydrogens = None

		# Calculated molecular weight
		self.weight = 0.0

		"""
		=== Processing Functions === 
		Functions for: Conversion to immutable types; Building of 
		neighbor tables; Calculation of Hybridization state, etc. 
		"""

		def immutable(table):
			"""Convert a list-based matrix or connection table into an
			immutable tuple of tuples."""
			if not table:
				return ()
			table = table[:]
			for i in range(len(table)):
				if type(table[i]) == str:
					continue
				table[i] = tuple(table[i])
			return tuple(table)

		def make_immutable(func):
			"""Function decorator version."""
			def wrap(arg1, arg2=None):
				if arg2:
					return immutable(func(arg1, arg2))
				return immutable(func(arg1))
			return wrap

		@make_immutable
		def generate_alpha_table(mat):
			"""Compute the alpha (direct neighbor) connection table 
			upfront from the connection matrix."""
			sz = len(mat)
			table = [[] for x in range(sz)]
			for i in range(sz):
				for j in range(sz):
					if i == j:
						continue
					if mat[i][j]:
						table[i].append(j)
			return table

		@make_immutable
		def generate_beta_table(alpha):
			"""Compute the beta (neighbor of neighbor) connection table
			upfront from the alpha table."""
			sz = len(alpha)
			table = [[] for x in range(sz)]
			for i in range(sz):
				for j in alpha[i]:
					for k in alpha[j]:
						if k == i:
							continue
						table[i].append(k)
			return table

		def compute_hybridizations(bondOrderMat, neighborTable):
			"""
			Generate Hybridization State for each atom.
			This is determined by analyzing the number of pi bond
			systems from the connection table. 

			Each atom is one of: {'sp', 'sp2', 'sp3', 'error'}
			"""
			HYBRID_VALUES = {0: 'sp3', 1: 'sp2', 2: 'sp'}

			sz = len(bondOrderMat)
			hybrids = ['error' for x in range(sz)]
			for i in range(sz):
				numPi = 0 # Number of pi systems
				for n in neighborTable[i]:
					bond = bondOrderMat[i][n]
					if bond >= 2:
						numPi += bond - 1
				if numPi in HYBRID_VALUES:
					hybrids[i] = HYBRID_VALUES[numPi]

			return tuple(hybrids)

		def compute_degrees(atomTypes, neighborTable):
			"""
			Compute the atom degrees. For carbon, this is the number of
			other carbons it is directly attached to. For other atoms,
			it is the degree of the carbon it is attached to.
			"""
			# FIXME: Degree of ethers, amides, esters, carbonyls?
			# Granted these are 'functional groups' and not lone atoms. 
			def carbon_degree(atom):
				deg = 0
				for n in neighborTable[atom]:
					if atomTypes[n].upper() == 'C':
						deg += 1
				return deg

			sz = len(atomTypes)
			degrees = [-1 for x in range(sz)]
			for i in range(sz):
				aType = atomTypes[i].upper()
				deg = 0
				if aType == 'C':
					deg = carbon_degree(i)
				else: 
					# Calculate degree for non-carbon atoms
					deg = -1
					for n in neighborTable[i]:
						# FIXME: Not sure what to do if two carbons, 
						# eg. the Oxygen in 'COC'. For now, just use 
						# the first carbon.
						if atomTypes[n].upper() == 'C':
							deg = carbon_degree(n) # FIXME: Redundant
							break
				degrees[i] = deg
			return tuple(degrees)

		def calculate_hydrogens(types, hybridizations, neighbors):
			"""
			Calculate number of hydrogens attached to each atom.
			Factors in hybridization state and number of neighbors.
			"""
			substituents = {'sp3': 4, 'sp2': 3, 'sp': 2}
			sz = len(types)
			hydrogens = [0 for x in range(sz)]

			for i in range(sz):
				# Number of subtituents is based on hybridization
				h = hybridizations[i]
				sub = 0
				if h in substituents:
					sub = substituents[h]

				# Minus number of current substituents
				sub -= len(neighbors[i])

				# Certain atoms prefer lone pair of electrons rather
				# than substituents (valency). 
				aType = types[i].upper()
				if aType == 'N':
					sub -= 1
				elif aType == 'O':
					sub -= 2
				elif aType in ['F', 'CL', 'BR']:
					sub -= 3

				if sub > 0:
					hydrogens[i] = sub
			return tuple(hydrogens)

		def calculate_molecular_weight(types, hydrogens):
			"""Calculate molecular weight."""
			weight = 0.0
			for atom in types:
				weight += ATOMIC_WEIGHTS[atom.upper()]

			hWeight = ATOMIC_WEIGHTS['H']
			for i in range(len(hydrogens)):
				weight += hydrogens[i] * hWeight
			return weight

		"""
		Setup from constructor input.
		"""

		# FIXME: Handle data that isn't supplied.

		self.connectMat = immutable(connectMat)
		self.bondOrderMat = immutable(bondOrderMat)

		self.alphaAtoms = generate_alpha_table(connectMat)
		self.betaAtoms = generate_beta_table(self.alphaAtoms)

		self.types = tuple(types)
		self.charges = tuple(charges)
		self.isotopes = tuple(isotopes)
		self.degrees = compute_degrees(types, self.alphaAtoms)
		self.hybridizations = compute_hybridizations(bondOrderMat, 
								self.alphaAtoms)

		self.hydrogens = calculate_hydrogens(types, self.hybridizations,
								self.alphaAtoms)
	
		self.weight = calculate_molecular_weight(types, self.hydrogens)

		self.smiles = smiles

	def setRings(self, rings):
		"""
		Set the rings that were found in ring perception.
		This sets certain flags.
		"""
		# Set ring membership flags and references
		self.isInRing = [False for x in range(self.size)]
		for ring in rings:
			for atom in ring:
				self.isInRing[atom] = True

		self.rings = rings

	def setRingGroups(self, ringGroups):
		"""
		Sets the ring groups that were found in ring perception.
		Sets certain flags an references.
		"""
		# Set ring group references
		def ring_group_membership(atom):
			for ringGroup in ringGroups:
				for ring in ringGroup:
					if atom in ring:
						return ringGroup
			return None

		self.ringGroupRef = [None for x in range(self.size)]
		for atom in range(self.size):
			self.ringGroupRef[atom] = ring_group_membership(atom)

	def setChains(self, chains):
		"""
		Set the chains that were found in chain perception.
		This sets certain flags and references.
		"""
		# Set chain membership flags and references
		# XXX: This does not include capping substituents!
		self.chainRef = [None for x in range(self.size)]
		self.isInChain = [False for x in range(self.size)]
		for chain in chains:
			for atom in chain:
				self.chainRef[atom] = chain
				self.isInChain[atom] = True

		self.chains = chains 

	def __setattr__(self, k, v):
		"""Limit the ability to manage the object's dictionary."""
		valid = ('size',
				'connectMat',
				'bondOrderMat',
				'types',
				'charges',
				'isotopes',
				'degrees',
				'hybridizations',
				'alphaAtoms',
				'betaAtoms',
				'ringSystem',
				'smiles',
				'informalName',
				'isInRing',
				'isInChain',
		)
		#if k not in valid:
		#	raise Exception("Cannot set as `%s`, invalid key." % k)
		### TODO/TEMP commented out
		###if k in self.__dict__ and self.__dict__[k]:
		###	raise Exception("Cannot reset `%s`." %k)

		self.__dict__[k] = v

	def __eq__(self, other):
		"""
		DANGEROUS! Compare two molecules to see if they are "equal".
		Right now, this is just a heuristic comparison and is not
		actual equality. Determining actual equality would require a 
		lot of work.
		"""
		if type(other) != Molecule:
			return False
		
		if self.size != other.size:
			return False

		if self.connectMat != other.connectMat:
			return False

		if self.bondOrderMat != other.bondOrderMat:
			return False

		# FIXME: This doesn't mean two molecules are the same. 
		# They could have different connection matrices for instance. 
		# They may also have the same layout, but different charges...
		# This is just a quick compare. (Harmful?)
		return True

	def __str__(self):
		"""String representation of the molecule for debugging."""
		txt = ""
		txt += "Name: %s\n" % str(self.informalName)
		txt += "Smiles: %s\n" % str(self.smiles)

		txt += "\nMolecular Weight: %f" % self.weight
		txt += "\nTypes: %s" % str(self.types)
		txt += "\nCharges: %s" % str(self.charges)
		txt += "\nIsotopes: %s" % str(self.isotopes)
		txt += "\nHybridizations: %s" % str(self.hybridizations)
		txt += "\nDegrees: %s" % str(self.degrees)
		txt += "\nHydrogens: %s" % str(self.hydrogens)

		l = [x for x in range(self.size) if self.isInRing[x]]
		txt += "\n\nRing atoms: %s" % l
		
		l = [x for x in range(self.size) if self.isInChain[x]]
		txt += "\nChain atoms: %s" % l

		txt += "\n\n"

		# XXX: Won't print >= 100 atoms nicely. Not that it would be
		# wise to print out such systems in the terminal...
		line1 = " "*5 if self.size < 10 else " "*6
		line2 = " "*5 if self.size < 10 else " "*6

		# Header atoms and header numbers
		for i in range(self.size):
			if len(self.types[i]) > 1:
				line1 += self.types[i]
			else:
				line1 += "%s " % self.types[i]
			if i < 10 or i %2 == 0:
				line2 += "%d " % i
			else:
				line2 += " "

		txt += "CONNECTION TABLE:"
		txt += "\n%s" % line1
		txt += "\n%s" % line2
		txt += "\n"

		def row_header(i):
			# FIXME: Can't print systems with > 99 atoms or print 
			# 3 char atoms... but need for either is very unlikely
			atom = self.types[i] 
			if len(atom) < 2:
				atom += " "
			if self.size < 10 or i >= 10:
				return "%s %d  " % (atom, i)
			return "%s %d   " % (atom, i) 

		# Bond Order Graph data
		for i in range(self.size):
			ln = row_header(i)
			for j in range(self.size):
				ln += str(int(self.bondOrderMat[i][j])) + " " \
						if self.bondOrderMat[i][j] else ". "
			txt += "%s\n" % ln

		# Lots of information 
		txt += "\nLabels; Hybridization; Degree; Alpha and Beta Atoms:\n"
		for i in range(self.size):
			ln = row_header(i)
			hybrid = str(self.hybridizations[i])
			hybrid = hybrid if hybrid != "error" else "err"
			hybrid += "  " if len(hybrid) == 3 else "   "
			degree = str(self.degrees[i])
			degree += "  " if len(degree) == 1 else " "
			ln += hybrid + degree

			ln += str(self.alphaAtoms[i])
			# FIXME: This only prints well in my Bash configuration AFAIK
			if len(ln) < 24:
				ln += "\t\t"
			else:
				ln += "\t"
			ln += str(self.betaAtoms[i])
			txt += "%s\n" % ln

		# DEBUG: Print second time 
		#txt += "\nInformal Name: %s" % str(self.informalName)
		#txt += "\nSmiles: %s" % str(self.smiles)

		return txt



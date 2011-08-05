"""
TODO: Doc.
"""

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
		"""

		# Number of atoms. Typically non-Hydrogen included.
		self.size = 0

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

		# Ring system
		self.ringSystem = None

		# Smiles text
		self.smiles = None

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

		"""
		Setup from constructor input.
		"""

		# FIXME: Handle data that isn't supplied.

		self.size = len(connectMat)

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

		self.smiles = smiles

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
				'smiles'
		)
		#if k not in valid:
		#	raise Exception("Cannot set as `%s`, invalid key." % k)
		if k in self.__dict__ and self.__dict__[k]:
			raise Exception("Cannot reset `%s`." %k)

		self.__dict__[k] = v

	def print_matrix(self):
		"""
		Print the molecular data for debugging.
		"""
		if self.smiles:
			print "===== Report for: %s =====\n" % self.smiles

		print "Types: %s" % str(self.types)
		print "Charges: %s" % str(self.charges)
		print "Isotopes: %s" % str(self.isotopes)
		print "Hybridizations: %s" % str(self.hybridizations)
		print "Degrees: %s" % str(self.degrees)
		print ""

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

		print line1
		print line2

		def row_header(i):
			# FIXME: Can't print systems with > 99 atoms or print 
			# 3 char atoms... but need for either is very unlikely
			atom = self.types[i] 
			if len(atom) < 2:
				atom += " "
			if self.size < 10 or i >= 10:
				return "%s %d  " % (atom, i)
			return "%s %d   " % (atom, i) 

		# Graph data
		for i in range(self.size):
			ln = row_header(i)
			for j in range(self.size):
				ln += str(int(self.bondOrderMat[i][j])) + " " \
						if self.bondOrderMat[i][j] else ". "
			print ln

		# Alpha atoms:
		print "\nAlpha Atoms:"
		for i in range(self.size):
			ln = row_header(i)
			ln += str(self.alphaAtoms[i])
			print ln

		# Beta atoms:
		print "\nBeta Atoms:"
		for i in range(self.size):
			ln = row_header(i)
			ln += str(self.betaAtoms[i])
			print ln

		if self.smiles:
			print "\n===== End report for: %s =====" % self.smiles


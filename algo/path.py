# Path Algorithms
#

class ShortestPaths(object):
	"""
	Floyd-Warshall shortest path between vertex pairs.

	Computes the shortest path between every vertex pair and stores the
	weights and path reconstruction matrix. This algorithm is extremely
	slow: O(n^3) worst case. Path reconstruction takes additional
	recursive computation between each vertex pair.

	You can do the calculation work upfront by calling calculate(), or
	wait for the first call of getWeight() or findPath() to do so.
	"""

	def __init__(self, matrix):
		"""
		Supply the graph adjacency matrix.
		Not connected must be represented by float('Infinity').
		"""
		self.matrix = matrix
		self.weights = None
		self.paths = None
		pass

	def size(self):
		"""
		Return the size of the matrix.
		"""
		return len(self.matrix)

	def calculate(self):
		"""
		Calculate the shortest paths.
		Call to perform calculations immediately, otherwise they will
		be calculated upon calling the other methods.
		"""
		# Do not calculate twice!
		if self.weights:
			return 		
		data = shortest_paths(self.matrix)
		self.weights = data[0]
		self.paths = data[1]

	def getWeight(self, i, j):
		"""
		Get the weight of the shortest path between i and j.
		If the result set has not already been calculated, it will now.
		"""
		# Floyd-Warshall. Results will be cached. 
		self.calculate()
		return self.weights[i][j]

	def findPath(self, i, j, include=False):
		"""
		Get the shortest path between i and j.
		If the result set has not already been calculated, it will now.
		Additionally, a recursive algorithm has to run in order to 
		reconstruct the path. 
		Returns one of:
			* A list containing the path
			* None, if (i->j) is an edge itself.
			* float('Infinity'), if DNE
		"""
		# Floyd-Warshall. Results will be cached. 
		self.calculate()

		# TODO: Mechanism to cache?
		ret = reconstruct_path(self.paths, i, j)

		if not include or type(ret) == float:
			return ret

		if include and ret == None:
			ret = []

		ret.insert(0, i)
		ret.append(j)
		return ret

def shortest_paths(matrix):
	"""
	Floyd-Warshall shortest paths algorithm.

	Constructs two matrices, one representing the cost between any two
	vertex pairs, the other including the path reconstruction matrix.
	Return value is a tuple including these.

	I should investigate using Johnson's Algorithm for sparse acyclic
	graphs.
	"""
	# Path does not exist for infinite values. 
	inf = float('Infinity')
	DNE = [inf, -inf]

	length = len(matrix)

	# Holds the length of the path between i and j as c[i][j]
	# We build this Dynamic Programming matrix as we increment k. 
	cost = [[inf for x in range(length)] for xx in range(length)]

	# Path reconstruction matrix.
	# For the shortest path between i and j, path[i][j] is the next
	# vertex along the path. Next take path[vert][j], and so forth...
	# (See the reconstruct path routine.)
	path = [[inf for x in range(length)] for xx in range(length)]

	# Initialize cost^(0) = matrix 
	for i in range(length):
		for j in range(length):
			cost[i][j] = matrix[i][j]
			if cost[i][j] not in DNE:
				path[i][j] = None # XXX: None denotes directly exists

	# Floyd's Algorithm
	for k in range(length):
		for i in range(length):
			for j in range(length):
				if i == j:
					continue
				if cost[i][j] > cost[i][k] + cost[k][j]:
					cost[i][j] = cost[i][k] + cost[k][j]
					path[i][j] = k

	return (cost, path)

def reconstruct_path(pathMat, i, j):
	"""
	Recursive algorithm to reconstruct the path between two vertices
	from the matrix that Floyd-Warshall returns.
	"""
	# FIXME: Don't keep redefining these recursively.
	# FIXME: Will recursion kill the Python stack?
	inf = float('Infinity')
	DNE = [inf, -inf]

	k = pathMat[i][j]

	if k in DNE:
		# Path does not exist between i and j. 
		# TODO: Throw exception?
		return inf
	if k == None:
		# This path is an edge itself. 
		return None

	l = reconstruct_path(pathMat, i, k)
	r = reconstruct_path(pathMat, k, j)

	cur = [k]
	if type(l) == list:
		cur = l + cur
	if type(r) == list:
		cur = cur + r

	return cur


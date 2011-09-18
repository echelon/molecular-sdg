
class Point(object):
	"""Represents a 2D position."""

	def __init__(self, x=None, y=None):
		self.x = x
		self.y = y

		if type(x) in [tuple, list]:
			self.x = x[0]
			self.y = x[1]

	def __eq__(self, o):
		"""Equal if within a certain delta."""
		DELTA = 0.00005
		if abs(self.x - o.x) > DELTA:
			return False
		if abs(self.y - o.y) > DELTA:
			return False
		return True

	def __str__(self):
		x = 'N' if self.x == None else str("%0.1f" % self.x)
		y = 'N' if self.y == None else str("%0.1f" % self.y)
		return "pt(%s, %s)" % (x, y)

	def __repr__(self):
		return str(self)


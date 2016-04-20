import math
import numpy as np
from munkres import Munkres

# define the Diagram class
# this is used to store persistence diagrams
class Diagram:
	def __init__(self):
		self.points = []
		self.infpoints = []
	def addpt(self, pt):
		self.points.append(pt)

	def addinfpt(self, pt):
		self.infpoints.append(pt)

	# add scale and translate methods later

def inf_pt_dist(diagram1, diagram2, q):
	"""
	This function computes the distance between the 
	points that correspond to essential classes,
	i.e. classes where the second component is infinite

	This function assumes that the sets of points are
	ordered by the first component (the finite one).
	This should be the case, as points are added to 
	the diagram in this order during construction.

	If the the diagrams have a different number of points
	at infinity, the distance between them is 
	infinite.

	Args:
		diagram1:  first persistence diagram
		diagram2:  second persistence diagram
		q:  we use L_q distance

	Return:  distance between points at infinity
		for the two diagrams

		"""

	n = len(diagram1.infpoints)
	m = len(diagram2.infpoints)

	if m != n:
		print("Warning:  different number of points are infinity")
	# this code needs to be optimized
	distance = 0
	for i in range(n):
		distance += math.pow(abs(diagram1.infpoints[i] - 
			diagram2.infpoints[i]), q)
	return distance

def L1dist(x,y):
	"""
	This function computes the L1 distance between
	two vectors x and y
	"""
	c = abs(a[0]-b[0]) + abs(a[1] - b[1])
	return c

def diag_len(x):
	"""
	This function calculates the distance from the point
	x (in the persistence diagram) to the diagonal.  This
	function uses the L1 distance and will be used in later 
	functions.
	"""
	return x[1] - x[0]


def finite_pt_dist(diagram1, diagram2, q):
	"""
	This function computes the smallest distance between the
	finite points of two persistence diagrams.  This requires
	matching each point in the first diagram to a point in the 
	second diagram or a point on the diagonal such that the distance
	between the pairs of points is minimized.

	The solution to this problem involves the Munkres (or Hungarian)
	algorithm, which has already been implemented in Python.

	This function returns the distanc between the two diagrams
	"""

	n = len(diagram1.points)
	m = len(diagram2.points)

	# if there are no points, the distance is zero
	if n + m == 0:
		return 0
	else:
		# this code can probabily be optimized
		dist_mat = []
		for i in range(n):
			row = []
			for j in range(m):
				length = L1dist(diagram1.points[i],diagram2.points[j])
				row.append(math.pow(length,q))
			a = math.pow(diag_len(diagram1.points[i]), q)
			for k in range(n):
				row.append(a)
			matrix.append(row)
		# create row with distance to diagonal
		row = []
		for j in range(m):
			row.append(math.pow(diag_len(diagram2.points[j]), q))
		for i in range(n):
			row.append(0)
		# add i copies of row to matrix
		for i in range(m):
			matrix.append(row)
		m = Munkres()
		# this extracts the indices of the shortest distance
		indices = m.compute(matrix)
		# now we can compute the total distance
		total_dist = 0
		for row, col in indices:
			value = matrix[row][column]
			total += value
		return total

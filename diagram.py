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


def read_mesh_graph(filename, d):
	file = open(filename, "r")
	line = file.readline()
	splitline = line.split()
	num_vert = int(splitline[0])
	num_edges = int(splitline[1])

	dict_vert = {}
	list_edges = []

	# dictionary of vertices {i: v_i}

	if d == 3:
		for i in range(num_vert):
			line = file.readline()
			splitline = line.split()

			dict_vert.update({i:(float(x) for x in splitline)})
	elif d == 2:
		for i in range(num_edges):
			line = file.readline()
			splitline = line.split()

			dict_vert.update({i:(float(x) for x in splitline)})

	for i in range(num_edges):
		line = file.readline()
		splitline = line.split()
		list_edges.append((int(splitline[0]),int(splitline[1])))

	list_vert = []
	for edge in list_edges:
		list_vert.append(e[0])
		list_vert.append(e[1])
	no_rep_vert = list(set(list_of_vert))
	dict_vert_clean = {}
	# does this need to be a set?
	# set in original code, list here
	for v in no_rep_vert:
		dict_vert_clean.update({v: dict_vert[v]})

	return(dict_vert_clean, list_edges)

def direction_order(dict_vert, list_edges, direction):
	"""
	This function computes the height of all vertices with
	respect to a specified direction and returns the height 
	of each vertex.

	The function also computes, for each vertex, a list of vertices
	that are lower (with respect to direction) 
	"""
	
	# this dictionary will contain the vertices that have
	# a lower height with respect to direction
	dict_neighbors = {}

	# this dictionary will map {v : height(v)} for
	# each vertex v
	list_heights = []
	dict_heights
	for v in dict_vert:
		coords = dict_vert[v]
		# this is height with respect to direction
		height = sum(vertex[i]*direction[i] for i in range(len(vertex)))
		dict_heights.update({v : height})
		list_heights.append((v,height))
	# sort vertices by height	
	list_heights.sort(key = lambda x: x[1])

	dict_neighbors = {v : [] for v in dict_vert.keys()}

	# for each vertex we find all other vertices that lower
	# height with respect to direction 
	for edge in list_edges:
		if dict_heights[edge[0]] > dict_heights[edge[1]]:
			# if first vertex is higher than second
			# add second vertex to list for first
			# vertex
			dict_neighbors[edge[0]].append(edge[1])

		elif dict_heights[edge[0]] == dict_heights[edge[1]]:
			# if they have the same height, we add the 
			# lower index to the list of the larger index
			if edge[0] > edge[1]:
				dict_heights[edge[0]].append(edge[1])
			# cannot have edge from vertex to itself
			# so if above condition is false, we must
			# have edge[1] > edge[0]
			else:
				dict_heights[edge[1]].append(edge[0])
		# if other conditions are not met, then edge[1] is
		# higher than edge[0]
		else:
			dict_neighbors[edge[1]].append(edge[0])
	return list_heights, dict_heights, dict_neighbors

# these classes are used for the Union-find algorithm
# this is a method we can use to detect connected class
# in the persistence diagrams

class Node:
	def __init__(self, index):
		self.parent = self
		self.index = index
	
	def changeroot(self, y):
		# this class is used to
		# update the root of tree of 
		# the node for x
		# this essentially merges two
		# trees
		self.parent = y

def Find(x):
	"""
	This function determines the root of the tree
	that x is in.  It works recursively.
	"""
	if x.parent == x:
		return x
	else:
		x.parent = Find(x.parent)
		return x.parent
def make_diagram(list_heights, dict_heights, dict_neighbors, 
		infinity = 100):
	"""
	This function makes a persistence diagrams given heights

	We use Union-Find algorithm to detect when there cycles,
	which represent merged classes
	"""
	# first we create empty an empty diagram class
	diagram = Diagram()
	
	dict_nodes = {}

	for i in list_heights:
		# here i is (v, height(v))
		v = i[0] # number of current vertex
		h = i[1] # height of current vertex
		if len(dict_neighbors[i[0]]) == 0:
			# if there are no points below the current point
			# it represents a new class and becomes
			# a note with itself as a parent
		else:
			list_compons = []
			for j in dict_neighbors[i[0]]:
				# Find(dict_nodes[j]).index returns the index
				# of the first node in the component in which
				# the vertex is lcoate
			set_compons = set(list_compons)
			ordered_compons = list(set_compons)
			
			if len(ordered_compons) == 1:
				node = Node(v)
				node.changeroot(dict_nodes[ordered_compons[0]])
				dict_nodes.update({v: node})
			else:
				list_birth_times = [dict_heights[j] for j in
						ordered_compons]
				birth = min(list_birth_times)
				count = 0

				ordered_compons.sort()
				for j in ordered_compons:
					if dict_heights[j] > birth:
						if dict_heights[j] < h - 1:
							diagram.addpt([dict_heights[j], h])
					if dict_heights[j] == birth:
						if count == 0:
							first_class = j
							node = Node(v)
							node.changeroot(dict_nodes[first_class])
							dict_nodes.update({v: node})
							count = 1
						else:
							if dict_heights[j] < h - 1:
								diagram.addpt([dict_heights[j],h])
				for k in ordered_compons:
					dict_nodes[k].changeroot(dict_nodes[first_class])
			final_compons = []
			for v in dict_nodes:
				final_compons.append(Find(dict_nodes[v]).index)
			set_final_compons = set(final_compons)
			for v in set_final_compons:
				diagram.addinfpt(dict_heights[v])

	return diagram



# these functions are used to sample directions
# the directions are used to construct persistence diagrams
# we need to construct persistence diagrams in many directions
# the 2-D functions will return exactly n directions, either
# randomly generated or evenly spaced
# the 3-d functions will return at least n directions
def sample_circle(n):
	"""
	This function gives evenly spaced directions
	on the cirlce and returns them in a 2 by n 
	matrix.  Each row represents a direction
	"""
	thetas = np.linspace(0, 2*np.pi, n)
	xs = np.cos(thetas)
	ys = np.sin(thetas)
	directions = np.array([xs,ys]).T
	return directions

def random_circle(n):
	"""
	This function randomly samples n points on the circle
	"""
	draws = np.random.uniform(low=0, high=2*np.pi, size=n)
	xs = np.cos(draws)
	ys = np.sin(draws)
	directions = np.array([xs,ys]).T
	return directions

def sample_sphere(n):
	"""
	This function returns approximately n points
	on the sphere.  It may return more if n is 
	not a perfect square.  It will always return at 
	least n points.  
	"""
	N = np.sqrt(n) + 1
	thetas = np.linspace(0, np.pi, N)
	phis = np.linspace(0, 2*np.pi, N)
	points = np.meshgrid(thetas, phis)
	xs = np.sin(points[0])*np.sin(points[1])
	ys = np.sin(points[0])*np.cos(points[1])
	zs = np.cos(points[0])
	
	matrix = np.vstack([xs.flatten(), ys.flatten(), zs.flatten()]).T
	return matrix


def random_sphere(n):
	"""
	This function returns approximately n randomly
	select points on the sphere.  It may return more if n is 
	not a perfect square.  It will always return at 
	least n points.  
	"""
	N = np.sqrt(n) + 1
	thetas = np.random.uniform(0, np.pi, N)
	phis = np.random.uniform(0, 2*np.pi, N)
	points = np.meshgrid(thetas, phis)
	xs = np.sin(points[0])*np.sin(points[1])
	ys = np.sin(points[0])*np.cos(points[1])
	zs = np.cos(points[0])
	
	matrix = np.vstack([xs.flatten(), ys.flatten(), zs.flatten()]).T
	return matrix

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

			dict_vert.update({i:list(float(x) for x in splitline)})
	elif d == 2:
		for i in range(num_edges):
			line = file.readline()
			splitline = line.split()
			dict_vert.update({i:list(float(x) for x in splitline)})
	for i in range(num_edges):
		line = file.readline()
		splitline = line.split()
		list_edges.append((int(splitline[0]),int(splitline[1])))
	list_vert  = []
	for edge in list_edges:
		list_vert.append(edge[0])
		list_vert.append(edge[1])
	no_rep_vert = list(set(list_vert))
	dict_vert_clean = {}
	# does this need to be a set?
	# set in original code, list here
	for v in no_rep_vert:
		dict_vert_clean.update({v: dict_vert[v-1]})
	file.close()
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
	dict_heights = {}
	for v in dict_vert.keys():
		coords = dict_vert[v]
		# this is height with respect to direction
		height = sum(coords[i]*direction[i] for i in range(len(coords)))
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
				dict_neighbors[edge[0]].append(edge[1])
			# cannot have edge from vertex to itself
			# so if above condition is false, we must
			# have edge[1] > edge[0]
			else:
				dict_neighbors[edge[1]].append(edge[0])
		# if other conditions are not met, then edge[1] is
		# higher than edge[0]
		else:
			dict_neighbors[edge[1]].append(edge[0])
	return list_heights, dict_heights, dict_neighbors

# these classes are used for the Union-find algorithm
# this is a method we can use to detect connected class
# in the persistence diagrams

class Tree:
	def __init__(self, name):
		self.parent = self
		self.name = name
		self.rank = 0
	
def Find(x):
	"""
	This function determines the root of the tree
	that x is in.  It works recursively. x should
	be an object of class Tree.
	"""
	if x.parent != x:
		x.parent = Find(x.parent)	
	return x.parent

def height_Union(x, y, dict_heights):
	"""
	This function takes the union of two nodes.
	It does this by changing the root one tree
	to be the root of the other tree.  It changes the
	root based on height.  The root becomes the node with
	the lowest height.  So the node that is born first
	becomes the root.  

	If the two roots have the same height, the lowest 
	number becomes the root.  For example, if we have vertex
	1 and vertix 3 at the same height, vertex 1 will become
	the root.  

	Inputs:
		x,y:  objects of Tree class
		dict_heights:  a map v-> h, where v is a vertex,
		which should be the .name of some tree and h is
		the height with respect to some direction
	"""
	x_root = Find(x)
	y_root = Find(y)
	if x_root == y_root:
		return None
	if dict_heights[x_root.name] < dict_heights[y_root.name]:
		y_root.parent = x_root
	elif dict_heights[x_root.name] == dict_heights[y_root.name]:
		if x.name < y.name:
			y_root.parent = x_root
		else:
			x_root.parent = y_root
	else:
		x_root.parent = y_root


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


# now we write a function that actually computes distances
def pairwise_distances(list_shapes, list_directions):
	"""
	This function actually calculates the pairwise distance between
	between all objects in list_shapes.  The distance is estimated
	using the vectors in list_directions.  It returns a pairwise
	distance matrix.  
	"""

	# first we calculate the persistence diagram for each object
	# in each direction
	l_diagrams = []
	for shape in list_shapes:
		shape_diagrams = []
		for direction in list_directions:
			l_heights, d_heights, d_n = direction_order(shape[0],
					shape[1], direction)
			shape_diagram = make_diagram(l_heights, d_heights,
					d_n)
			shape_diagrams.append(shape_diagram)
		l_diagrams.append(shape_diagrams)
	N = len(list_shapes)
	K = len(list_directions)
	dists = np.empty((N,N))
	for i in range(N):
		for j in range(i,N):
			finite_dist = 0
			infinite_dist = 0
			for k in range(K):
				finite_dist += finite_pt_dist(l_diagrams[i][k],
						l_diagrams[j][k], 1)
				infinite_dist += inf_pt_dist(l_diagrams[i][k],
						l_diagrams[j][k], 1)
			dists[i,j] += finite_dist + infinite_dist
	dists += dists.T
	return dists

def make_diagram(dict_heights, dict_neighbors):
	dict_trees = {v: Tree(v) for v in dict_heights.keys()}
	vertices = set(dict_heights.keys())
	heights = set(dict_heights.values())
	heights_sorted = sorted(list(heights))
	dict_below_h = {h: [] for h in heights}
	for h in heights:
		for v in vertices:
			if dict_heights[v] <= h:
				dict_below_h[h].append(v)
	dict_deaths = {v: None for v in dict_heights.keys()}
	dead = set()
	for v in vertices:
		if len(dict_neighbors[v]) != 0:
			dict_deaths.update({v: dict_heights[v]})
			dead.add(v)
	alive = vertices.difference(dead)
	for h in heights_sorted:
		died = set()
		for v in dict_below_h[h]:
			for neigh in dict_neighbors[v]:
				height_Union(dict_trees[v], dict_trees[neigh],
						dict_heights)
		for v in alive:
			if dict_neights[v] <= h:
				own_root = v != Find(dict_trees[v]).name
				if own_root:
					died.add(v)
					dict_deaths.update({v: h})
		alive = alive.difference(died)
	diag = diagram.Diagram()
	for v in dict_heights.keys():
		if dict_deaths[v] == None:
			diag.addinfpt(dict_heights[v])
		else:
			diag.addpt((dict_heights[v], dict_deaths[v]))
	
	return diag

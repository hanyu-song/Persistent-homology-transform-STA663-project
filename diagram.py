import math
import multiprocessing as mp
import numpy as np
import scipy.io as sio
import glob
from scipy.optimize import linear_sum_assignment
from munkres import Munkres
from numba import jit

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
	c = abs(x[0]-y[0]) + abs(x[1] - y[1])
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
		dist_mat = make_dist_mat(diagram1, diagram2, q)
		# now we can compute the total distance
		row_ind, col_ind = linear_sum_assignment(dist_mat)
		total_dist = 0
		for i in range(len(row_ind)):
			row = row_ind[i]
			col = col_ind[i]
			value = dist_mat[row][col]
			total_dist += value
		return total_dist


def read_mesh_graph(filename):
	file = open(filename, "r")
	line = file.readline()
	splitline = line.split()
	num_vert = int(splitline[0])
	num_edges = int(splitline[1])

	dict_vert = {}
	list_edges = []

	# dictionary of vertices {i: v_i}

	for i in range(num_vert):
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
			shape_diagram = make_diagram(d_heights, d_n)
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
			dists[i,j] += (finite_dist + infinite_dist)/K
	dists += dists.T
	return dists

def scaled_distance(list_objects, matrix_dir):
	# number of objects
	N = len(list_objects)
	# number of directions	
	k = matrix_dir.shape[1]
	# compute centering constant
	K = 0
	for i in range(k):
		K += np.cos(2*np.pi*i/k)**2
	for i in range(N):
		obj = list_objects[i][0]
		prod = obj.dot(matrix_dir)
		# minimum of eac column is lambda_i
		lambdas = prod.min(axis = 0)
		u = (1/K)*matrix_dir.dot(lambdas)
		list_objects[i][0] = obj -u[np.newaxis,:]
		# now we scale each object
		L = -lambdas.sum()
		list_objects[i][0] /= L
	# now we create a list of diagrams
	# each item in the list is a list of
	# diagrams for one object in each direction
	# so it is a list of length N with sublists
	# of size k
	l_diagrams = []
	for i in range(N):
		shape = list_objects[i]
		# this is will hold all diagrams for
		# the current object
		shape_diagrams = []
		for j in range(k):
			direction = matrix_dir[:,j].T
			# now we transform the data into forms we
			# can use with our existing functions
			# to do this, we make a dictionary mapping
			# each vertex to its coordinates
			num_rows = shape[0].shape[0]
			num_edges = shape[1].shape[0]
			dict_vert = {i+1: shape[0][i,:] for i in range(num_rows)}
			list_edges = [list(shape[1][i,:]) for i in range(num_edges)]
			# make the diagram for jth direction
			l_heights, d_heights, d_n = direction_order(dict_vert,
					list_edges, direction)
			shape_diagram = make_diagram(d_heights, d_n)
			shape_diagrams.append(shape_diagram)
		l_diagrams.append(shape_diagrams)
	dist_mat = np.zeros((N,N))
	print("Beginning to calculate distances")
	for i in range(N):
		for j in range(i+1, N):
			print("Calculating distance between", i, "and", j)
			list_dists = []
			# to find the actual distance between two objects
			# we need to consider the distance under various
			# rotations.  rotations, however, are simply other
			# persistence diagrams.  Thus we calculate the distance
			# between a diagram and all other persistence diagrams
			# for the second diagram, which can be seen as 
			# rotations of the second diagram
			for shift in range(k):
				#print("Considering shift: ", shift)
				finite_dist = 0
				infinite_dist = 0
				for cur in range(k):
					#	print("Calculating finite distance for", cur)
					finite_dist += finite_pt_dist(l_diagrams[i][cur],
							l_diagrams[j][(cur + shift)% k],
							1)
					infinite_dist += inf_pt_dist(l_diagrams[i][cur],
							l_diagrams[j][(cur + shift)% k],
							1)
				list_dists.append(finite_dist + infinite_dist)
			actual_dist = min(list_dists)
			dist_mat[i,j] += actual_dist/k
			print(actual_dist/k)
			print(calculate_distance(i, j, k, l_diagrams))
	dist_mat += dist_mat.T
	return dist_mat



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
			if dict_heights[v] <= h:
				own_root = v != Find(dict_trees[v]).name
				if own_root:
					died.add(v)
					dict_deaths.update({v: h})
		alive = alive.difference(died)
	diag = Diagram()
	for v in dict_heights.keys():
		if dict_deaths[v] == None:
			diag.addinfpt(dict_heights[v])
		else:
			# only add point if not on the diagonal
			# this is designed to reduce the complexity
			# of the problem for the Munkres algorithm
			if dict_heights[v] != dict_deaths[v]:
				diag.addpt((dict_heights[v], dict_deaths[v]))
	
	return diag




# break up distance functions
def make_dist_mat(diagram1, diagram2, q):
	n = len(diagram1.points)
	m = len(diagram2.points)
	dist_mat = np.zeros((n+m,n+m))
	for i in range(n):
		for j in range(m):
			# distance from point i in diagram 1 to joint j in diagram 2

			distance = pow(L1dist(diagram1.points[i], diagram2.points[j]),q)
			dist_mat[i, j] += distance
		# distance from point i of diagram 1 to diagonal
		dist_to_diag = pow(diag_len(diagram1.points[i]), q)
		dist_mat[i, m:][np.newaxis, :] += dist_to_diag*np.ones((1,n))
		# now we consider the m copies of the diagonal
	for j in range(m):
		# distiance from jiagonal to point j of diagram 2
		dist_to_diag = pow(diag_len(diagram2.points[j]),q)
		dist_mat[n:,j][:,np.newaxis] += dist_to_diag*np.ones((m,1))
	return dist_mat

@jit
def make_dist_jit(diagram1, diagram2, q):
	n = len(diagram1.points)
	m = len(diagram2.points)
	dist_mat = np.zeros((n+m,n+m))
	for i in range(n):
		for j in range(m):
			# distance from point i in diagram 1 to joint j in diagram 2

			distance = pow(L1dist(diagram1.points[i], diagram2.points[j]),q)
			dist_mat[i, j] += distance
		# distance from point i of diagram 1 to diagonal
		dist_to_diag = pow(diag_len(diagram1.points[i]), q)
		for j in range(m):
			dist_mat[i, j] += dist_to_diag
		# now we consider the m copies of the diagonal
	for j in range(m):
		# distiance from jiagonal to point j of diagram 2
		dist_to_diag = pow(diag_len(diagram2.points[j]),q)
		dist_mat[n:,j][:,np.newaxis] += dist_to_diag*np.ones((m,1))
	return dist_mat


def calculate_distance(i, j, k, l_diagrams):
	"""
	This function calculates the distance between i and j modulo rotation.
	To do this, we calculate the distance between i and all rotations of j
	and take the minimal distance.  We do not actually rotate the object.
	Instead, we realize that the persistence diagram of a rotation is
	is persistence diagram of a different direction.  Thus we simply rename the 
	persistence diagrams.  

	Args:
		i:  index of first object
		j:  index of second object
		k:  number of directions
		l_diagrams:  a list of diagrams for each object in each of the
			k directions
	
	Return:
		The returns a tuple with 3 number: i, j, and the distance between
		these objects under rotation.  The return of i and j is for convenience
		in parallelization.
	"""
	list_dists = []
	for shift in range(k):
		finite_dist = 0
		infinite_dist = 0
		for cur in range(k):
			finite_dist += finite_pt_dist(l_diagrams[i][cur],
										  l_diagrams[j][(cur + shift)%k],
										  1)
			infinite_dist += inf_pt_dist(l_diagrams[i][cur],
										 l_diagrams[j][(cur + shift)%k],
										 1)
		list_dists.append(finite_dist + infinite_dist)
	actual_dist = min(list_dists)/k
	return i, j, actual_dist




def scaled_distance_par(list_objects, matrix_dir, workers):
	# number of objects
	N = len(list_objects)
	# number of directions	
	k = matrix_dir.shape[1]
	# compute centering constant
	K = 0
	for i in range(k):
		K += np.cos(2*np.pi*i/k)**2
	for i in range(N):
		obj = list_objects[i][0]
		prod = obj.dot(matrix_dir)
		# minimum of eac column is lambda_i
		lambdas = prod.min(axis = 0)
		u = (1/K)*matrix_dir.dot(lambdas)
		list_objects[i][0] = obj -u[np.newaxis,:]
		# now we scale each object
		L = -lambdas.sum()
		list_objects[i][0] /= L
	# now we create a list of diagrams
	# each item in the list is a list of
	# diagrams for one object in each direction
	# so it is a list of length N with sublists
	# of size k
	l_diagrams = []
	for i in range(N):
		shape = list_objects[i]
		# this is will hold all diagrams for
		# the current object
		shape_diagrams = []
		for j in range(k):
			direction = matrix_dir[:,j].T
			# now we transform the data into forms we
			# can use with our existing functions
			# to do this, we make a dictionary mapping
			# each vertex to its coordinates
			num_rows = shape[0].shape[0]
			num_edges = shape[1].shape[0]
			dict_vert = {i+1: shape[0][i,:] for i in range(num_rows)}
			list_edges = [list(shape[1][i,:]) for i in range(num_edges)]
			# make the diagram for jth direction
			l_heights, d_heights, d_n = direction_order(dict_vert,
					list_edges, direction)
			shape_diagram = make_diagram(d_heights, d_n)
			shape_diagrams.append(shape_diagram)
		l_diagrams.append(shape_diagrams)
	dist_mat = np.zeros((N,N))
	print("Beginning to calculate distances")
	list_inputs = [(i, j, k, l_diagrams) for i in range(N)
										 for j in range(i+1, N)]
	with mp.Pool(processes=workers) as pool:
		distances = pool.starmap(calculate_distance, list_inputs)
	for tup in distances:
		dist_mat[tup[0],tup[1]] = tup[2]
	dist_mat += dist_mat.T
	return dist_mat



def scaled_distance_short(list_objects, k, workers):
	# number of objects
	N = len(list_objects)
	# number of directions	
	matrix_dir = sample_circle(k).T
	# compute centering constant
	scaled_objects = center_scale(list_objects, matrix_dir)
	# now we create a list of diagrams
	# each item in the list is a list of
	# diagrams for one object in each direction
	# so it is a list of length N with sublists
	# of size k
	l_diagrams = make_diagrams(scaled_objects, matrix_dir)
	dist_mat = np.zeros((N,N))
	print("Beginning to calculate distances")
	list_inputs = [(i, j, k, l_diagrams) for i in range(N)
										 for j in range(i+1, N)]
	with mp.Pool(processes=workers) as pool:
		distances = pool.starmap(calculate_distance, list_inputs)
	for tup in distances:
		dist_mat[tup[0],tup[1]] = tup[2]
	dist_mat += dist_mat.T
	return dist_mat

def center_scale(list_objects, matrix_dir):
	k = matrix_dir.shape[1]
	# calculate scaling constant
	N = len(list_objects)
	K = 0
	for i in range(k):
		K += np.cos(2*np.pi*i/k)**2
	for i in range(N):
		obj = list_objects[i][0]
		prod = obj.dot(matrix_dir)
		# minimum of each column is lambda_i
		lambdas = prod.min(axis = 0)
		u = (1/K)*matrix_dir.dot(lambdas)
		list_objects[i][0] = obj - u[np.newaxis, :]
		# now we scale each object
		L = -lambdas.sum()
		list_objects[i][0] /= L
	return list_objects

def make_diagrams(list_objects, matrix_dir):
	k = matrix_dir.shape[1]
	N = len(list_objects)
	l_diagrams = []
	for i in range(N):
		shape = list_objects[i]
		shape_diagrams = []
		for j in range(k):
			direction = matrix_dir[:,j].T
			# now we transform the data into forms we
			# can use with our existing functions
			# to do this, we make a dictionary mapping
			# each vertex to its coordinates
			num_rows = shape[0].shape[0]
			num_edges = shape[1].shape[0]
			dict_vert = {i+1: shape[0][i,:] for i in range(num_rows)}
			list_edges = [list(shape[1][i,:]) for i in range(num_edges)]
			# make the diagram for jth direction
			l_heights, d_heights, d_n = direction_order(dict_vert,
					list_edges, direction)
			shape_diagram = make_diagram(d_heights, d_n)
			shape_diagrams.append(shape_diagram)
		l_diagrams.append(shape_diagrams)
	return l_diagrams

def read_closed_shapes(directory):
	"""
	This function reads in all .mat files a specified directory
	"""
	query = directory + "*.mat"
	files = glob.glob(query)
	shapes = []
	for file in files:
		vertices = sio.loadmat(file)['x']
		N = vertices.shape[0]
		edges = np.zeros((N,2))
		for i in range(N-1):
			edges[i,:] = np.array([i+1, i+2])
		shapes.append([vertices, edges])
	return shapes

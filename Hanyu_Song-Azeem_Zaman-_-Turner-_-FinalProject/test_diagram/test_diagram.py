import diagram
import pytest
import numpy as np

def test_on_diag():
	point = [1,1]
	dist_to_diag = diagram.diag_len(point)
	assert dist_to_diag == 0

def test_dist_to_diag():
	point = [1,2]
	dist_to_diag = diagram.diag_len(point)
	assert dist_to_diag == 1

def test_L1_dist():
	point1 = [2,1]
	point2 = [1,1]
	l1_dis = diagram.L1dist(point1,point2)
	assert l1_dis == 1

def test_dist_mat():
	diag1 = diagram.Diagram()
	diag1.addpt((0,0))
	diag2 = diagram.Diagram()
	diag2.addpt((1,1))
	dist_mat = diagram.make_dist_mat(diag1, diag2, 1)
	ans = np.array([[2,0],[0,0]])
	assert np.array_equal(dist_mat, ans)
	
def test_diag():
	objects = diagram.read_closed_shapes("../small_test/")
	dist_mat = diagram.distance_unscaled(objects, 10)
	N = len(objects)
	sum_diag = 0
	for i in range(N):
		sum_diag += dist_mat[i,i]
	assert sum_diag == 0

def test_sym():
	objects = diagram.read_closed_shapes("../small_test/")
	dist_mat = diagram.distance_unscaled(objects, 10)
	N = len(objects)
	assert (dist_mat == dist_mat.T).all()

def test_mc_unscaled():
	objects = diagram.read_closed_shapes("../small_test/")
	dist_mat = diagram.distance_unscaled(objects, 10)
	dist_mat_mc = diagram.distance_unscaled_mc(objects, 10, 4)
	assert np.array_equal(dist_mat, dist_mat_mc)

def test_mc_scaled():
	objects = diagram.read_closed_shapes("../small_test/")
	dist_mat = diagram.distance_scaled(objects, 10)
	objects = diagram.read_closed_shapes("../small_test/")
	dist_mat_mc = diagram.distance_scaled_mc(objects, 10, 4)
	assert np.array_equal(dist_mat, dist_mat_mc)

def test_dist_to_self_unscaled():
	objects = diagram.read_closed_shapes("../small_test/self_dist/")
	dist_mat = diagram.distance_unscaled(objects, 10)
	zeros = np.zeros((len(objects),len(objects)))
	assert np.array_equal(dist_mat, zeros)

def test_dist_to_self_scaled():
	objects = diagram.read_closed_shapes("../small_test/self_dist/")
	dist_mat = diagram.distance_scaled(objects, 10)
	zeros = np.zeros((len(objects),len(objects)))
	assert np.array_equal(dist_mat, zeros)

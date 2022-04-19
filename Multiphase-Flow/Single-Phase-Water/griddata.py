import numpy as np
from numpy import meshgrid
import math as m
from numpy.linalg import solve

def biharmonic_spline_interpolate(x,y,f,X,Y):
	#function of interpolation data
	#same with matlab biharmonic spline interpolate v4
	Z = np.zeros((X.shape[0], X.shape[1]))
	length = len(f)
	G = np.zeros((length, length), dtype = 'float32')
	for i in range(length):
		for j in range(length):
			if i != j:
				Ma = m.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2)
				if Ma >= 1e-7:
					G[i,j] = Ma**2*(np.log(Ma) -1)

	a = np.asarray(solve(G,f)).reshape(-1,1)
	g = np.zeros((a.shape[0], a.shape[1]))
	for i in range(Z.shape[0]):
		for j in range(Z.shape[1]):
			for k in range(length):
				Ma1 = m.sqrt((X[i,j] - x[k])**2 + (Y[i,j] -y[k])**2)
				if Ma1 >= 1e-7:
					g[k] = Ma1**2*(np.log(Ma1) -1)
				else:
					g[k] = Ma1**2 * (-100)
			Z[i,j] = sum(np.multiply(g,a))

	return Z



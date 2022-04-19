import numpy as np

def boundaries(X):
	'''This function will make zero values outside a matrix'''
	row = X.shape[0]
	col = X.shape[1]
	X_new = np.zeros((row+2, col+2))
	for i in range(row):
		for j in range(col):
			X_new[i+1, j+1] = X[i,j]
	return X_new



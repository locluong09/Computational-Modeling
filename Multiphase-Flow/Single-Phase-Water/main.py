import numpy as np
from griddata import biharmonic_spline_interpolate
import matplotlib.pyplot as plt
import math as m
from boundaries_zero import boundaries


#load data from txt file
x_coordinate = np.loadtxt(open('x.txt'))
y_coordinate = np.loadtxt(open('y.txt'))

thickness = np.loadtxt(open('thickness.txt'))
top = np.loadtxt(open('top.txt'))
Perm = np.genfromtxt(open('perm1.txt'),delimiter=",", dtype = np.float32)

#griddata for interpolate the perm, top and thickness
[X,Y] = np.meshgrid(x_coordinate, y_coordinate)
Thickness = biharmonic_spline_interpolate(thickness[:,0], thickness[:,1], 
											thickness[:,2], X, Y)
Top = biharmonic_spline_interpolate(top[:,0], top[:,1], top[:,2], X, Y)



#assert values outside the boundaries is zero
for i in range(Thickness.shape[0]):
	for j in range(Thickness.shape[1]):
		if Thickness[i,j] < 0:
			Thickness[i,j] = 0
		if Perm[i,j] < 0:
			Perm[i,j] = 0

#assert values of top outside the boundaries is 3100
for i in range(Top.shape[0]):
	for j in range(Top.shape[1]):
		if Top[i,j] > 3100:
			Top[i,j] = 3100
		#
		Top[0:3,0] = 3100
		Top[0:2,0:2] = 3100
		#
		Top[0:3, 26:28] = 3100
		Top[3,27] = 3100
		#
		Top[12:14, 26:28] = 3100
		Top[10:12, 26:28] = 3100
		#
		Top[0,0:28] = 3100
		Top[13,0:28] = 3100
		Top[0:14,0] = 3100
		Top[0:14,27] = 3100

''' I gonna use block center method to define each block. So I will define new_perm, new_thickness
and new_top by averaging 4 corner point values'''
Perm_new = np.zeros((Perm.shape[0]-1, Perm.shape[1] - 1))
Thickness_new = np.zeros((Thickness.shape[0]-1, Thickness.shape[1] -1))
Top_new = np.zeros((Top.shape[0]-1, Top.shape[1] -1))
for i in range(Perm.shape[0] -1):
	for j in range(Perm.shape[1]-1):
		Perm_new[i,j] = (Perm[i,j] + Perm[i+1,j] + Perm[i,j+1] + Perm[i+1, j+1])/4
		Thickness_new[i,j] = (Thickness[i,j] + Thickness[i+1,j] + Thickness[i,j+1] + Thickness[i+1, j+1])/4
		Top_new[i,j] = (Top[i,j] + Top[i+1,j] + Top[i,j+1] + Top[i+1, j+1])/4

#Create matrix with zeros boundaries outside
k = boundaries(Perm_new)
h = boundaries(Thickness_new)
Top = boundaries(Top_new)

#Create dx and dy matrix
deltax = np.zeros((Perm_new.shape[0], Perm_new.shape[1]))
deltay = np.zeros((Perm_new.shape[0], Perm_new.shape[1]))
x = x_coordinate
dx = np.zeros((1,len(x)-1))
for i in range(len(x)-1):
	dx[0,i] = x[i+1] - x[i]
y = y_coordinate
dy = np.zeros((len(y)-1,))
for i in range(len(y)-1):
	dy[i,] = y[i+1] - y[i]

for i in range(deltax.shape[0]):
	deltax[i,:] = dx
for j in range(deltay.shape[1]):
	deltay[:,j] = dy
#Zeros boundaties to dx dy
dx = boundaries(deltax)
dy = boundaries(deltay)

#~~~Order the active block by permeability~~~
Row,Col = k.shape[0],k.shape[1]
count = 0
order = np.zeros((Row,Col))
for i in range(Row):
	for j in range(Col):
		if k[i,j] > 0:
			count += 1
			order[i,j] = count
min_value = 1
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			if (h[i,j] < min_value) & (h[i,j] > 0):
				min_value = h[i,j]

for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			if h[i,j] == 0.0:
				h[i,j] = min_value

cou = 0
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			if k[i,j] == 0.0:
				cou += 1

cou1 = 0
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			if h[i,j] == 0.0:
				cou1 += 1

#Fluid properties
rho = 62.4 #lb/ft^3
B = 1 #RB/STB
miu = 1 #cP
beta_c = 1.127e-3

#Cross-sectional area
Ax = np.multiply(dy,h)
Ay = np.multiply(dx,h)


#Initialize N,S,W,E matrix 
N = np.zeros((Row,Col))
S = np.zeros((Row,Col))
W = np.zeros((Row,Col))
E = np.zeros((Row,Col))

#Define N,S,W,E and C matrix of transmisibility
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			N[i,j]=2*beta_c/(miu*B*(dy[i-1,j]/(Ay[i-1,j]*k[i-1,j])+ dy[i,j]/(Ay[i,j]*k[i,j])))
			S[i,j]=2*beta_c/(miu*B*(dy[i+1,j]/(Ay[i+1,j]*k[i+1,j])+ dy[i,j]/(Ay[i,j]*k[i,j])))
			W[i,j]=2*beta_c/(miu*B*(dx[i,j-1]/(Ax[i,j-1]*k[i,j-1])+ dx[i,j]/(Ax[i,j]*k[i,j])))
			E[i,j]=2*beta_c/(miu*B*(dx[i,j+1]/(Ax[i,j+1]*k[i,j+1])+ dx[i,j]/(Ax[i,j]*k[i,j])))

#change nan in matrix to zero
N[np.isnan(N)] = 0
S[np.isnan(S)] = 0
W[np.isnan(W)] = 0
E[np.isnan(E)] = 0

C = -(N+S+W+E)


#define Right hand side of the equation
Q = np.zeros((count,1))
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			Q[int(order[i,j])-1] = N[i,j]*Top[i-1,j]+S[i,j]*Top[i+1,j]+ W[i,j]*Top[i,j-1]+E[i,j]*Top[i,j+1]+C[i,j]*Top[i,j] 
Q = Q*rho/144

#Well block co-ordinate and properties
well_coor = [[4,6], [9,6], [4,20], [8,10], [12,13], [6,14]]
well_coor = np.asarray(well_coor)
well_spec = [[2, 1500], [2, 1500], [1, 14.7], [1, 3000], [1, 14.7], [1, 14.7]] #1-psf 2-q
well_ppt = [[0.25, -1], [0.25, -1], [0.25, -1], [0.25, -1], [0.25, -1], [0.25, -1]]
well_spec = np.asarray(well_spec)
well_ppt = np.asarray(well_ppt)

#Initialize omega
omega = np.zeros((Row,Col))
for i in range(well_coor.shape[0]):
	ii = well_coor[i,0] #index i of each well
	jj = well_coor[i,1] #index j of each well
	rw = well_ppt[i,0] #well_bore radius
	s = well_ppt[i,1] #skin
	re = 0.14*m.sqrt(dx[ii,jj]**2 + dy[ii,jj]**2)
	omega[ii,jj] = 2*m.pi*beta_c*k[ii,jj]*h[ii,jj]/(miu*B*(np.log(re/rw) +s ))

	if int(well_spec[i,0]) == 1:
		C[ii,jj] = C[ii,jj] - omega[ii,jj]
		Q[int(order[ii,jj])-1] = Q[int(order[ii,jj])-1] - omega[ii,jj]*well_spec[i,1] #for RHS
	else:
		Q[int(order[ii,jj])-1] = Q[int(order[ii,jj])-1] -well_spec[i,1]


'''Construct coefficient matrix of LHS for solving pressure distribution'''
LHS = np.zeros((count,count))
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			LHS[int(order[i,j])-1, int(order[i,j])-1] = C[i,j]
			if order[i-1,j] > 0:
				LHS[int(order[i,j])-1, int(order[i-1,j])-1] = N[i,j]
			if order[i+1,j] > 0:
				LHS[int(order[i,j])-1, int(order[i+1,j])-1] = S[i,j]
			if order[i,j-1] > 0:
				LHS[int(order[i,j])-1, int(order[i,j-1])-1] = W[i,j]
			if order[i,j+1] > 0:
				LHS[int(order[i,j])-1, int(order[i,j+1])-1] = E[i,j]

from numpy.linalg import lstsq, solve
RHS = Q

print('Number of block in reservoir:', (Row-2)*(Col-2))
print('Number of active block', count)

x = solve(LHS,RHS)
p_list = []
for i in range(well_coor.shape[0]):
	ii = int(well_coor[i,0])
	jj = int(well_coor[i,1])
	p_list.append(order[ii,jj])
	print("Well block number",order[ii,jj])
	#print(C[ii,jj])
	print("Well productivity index is:",omega[ii,jj])

p_list = map(int, p_list)
#print(list(p_list))
print("Well-block pressure:")
for item in p_list:

	print(x[item - 1])
	#print(RHS[item - 1])

Q_check = [i for i in np.arange(6)]
for i in range(well_coor.shape[0]):
	ii = well_coor[i,0]
	jj = well_coor[i,1]
	if int(well_spec[i,0]) ==1:
		Q_check[i] = -omega[ii,jj]*(x[int(order[ii,jj]) - 1] - well_spec[i,1])[0]
	else:
		Q_check[i] = well_spec[i,1]

print("Flow rates vector:",Q_check)
print("Material balance check",sum(Q_check))
print("Residual check at each well block")
print("Well 1",np.reshape(LHS[47,:],(1,-1)).dot(x) - Q[47,0])
print("Well 2",np.reshape(LHS[161,:],(1,-1)).dot(x) - Q[161,0])
print("Well 3",np.reshape(LHS[61,:],(1,-1)).dot(x) - Q[61,0])
print("Well 4",np.reshape(LHS[158,:],(1,-1)).dot(x) - Q[158,0])
print("Well 5",np.reshape(LHS[258,:],(1,-1)).dot(x) - Q[258,0])
print("Well 6",np.reshape(LHS[108,:],(1,-1)).dot(x) - Q[108,0])

print("Residual check all block:", LHS*x - RHS)
Pre_dis = np.zeros((Row,Col))
for i in range(Row):
	for j in range(Col):
		if order[i,j] > 0:
			Pre_dis[i,j] = x[int(order[i,j])-1,0]


plt.imshow(Pre_dis, cmap ='plasma')
plt.show()


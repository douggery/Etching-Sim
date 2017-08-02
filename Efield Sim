import matplotlib.pyplot as plt
import numpy as np
import math
import time

start=time.clock()

h=1; #distance in real space between grid points
l=60;
w=60; #number of grid points on each side of the rectangular mesh
error_limit=10**-6
epsilon_zero=1 #8.854*10**-8;

#calculate w, the residual weighting
t=math.cos(math.pi/l)+math.cos(math.pi/w);
omega=(8-math.sqrt(64-16*t**2))/(t**2);

#create empty boundary condition mesh
BC=np.zeros((l,w));

def generate_dirichlet_boundaries():
	#initialize boundary condition mesh where "1" is 
	#immutability and "0" is mutability

	#make boundaries have zero potential

	for i in range(l):
		BC[i][0]=1
		BC[i][w-1]=1
	
	for i in range(w):
		BC[0][i]=1
		BC[l-1][i]=1

	#generate boundaries from initial voltage and charge mesh
	
	for i in range(l):
		for j in range(w):
			if V[i][j] != 0:
				BC[i][j]=1
	print(BC)
	return BC

#initialize voltage mesh where the initial mesh 
#defines the boundary condition values and all else is void
V=np.zeros((l,w));
#for i in range(w//2):
#	V[l//2+10][i+w//4]=1
#	V[l//2-10][i+w//4]=-1

def create_voltage_plate(x_start,x_finish,y_start,y_finish,volts):
	#for simplicity, i take the origin to be the top left
	#of the image ; in principle it can be modified to have 
	#the origin in the bottom left with means the values [l-y][x]
	lx=(x_finish-x_start)
	ly=(y_finish-y_start)
	ls=math.sqrt(lx**2+ly**2)
	slope=ly/(0.000001+lx)
	if lx>ly:
		for i in range(lx):
			V[y_start+(i*slope)//1][i+x_start]=volts
	else:
		for j in range(ly):
			V[y_start+j][x_start+(j/slope)//1]=volts
	
create_voltage_plate(w//4,3*w//4,l//2-5,l//2-5,1)
create_voltage_plate(w//4,3*w//4,l//2+5,l//2+5,-1)
#create_voltage_plate(75,70,89,97,-5)
#create_voltage_plate(44,68,40,60,-5)
	
	
	
	
#V[0][0]=1
#V[0][9]=-1

#V[1][0]=1
#V[1][9]=-1

#V[2][0]=1
#V[2][9]=-1
	
print(V)

#initialize charge mesh
rho=np.zeros((l,w));

def create_point_charge(x_pos,y_pos,charge):
	rho[y_pos][x_pos]=charge
	
#create_point_charge(50,50,30)
#create_point_charge(80,80,-30)

#initialize residual mesh
R=np.zeros((l,w));

#calculate residual mesh
def calc_bulk_residuals():
	for i in range(1,l-1):
		for j in range(1,w-1):
			if BC[i][j]==1:
				R[i][j]=0
			else:
				a0=D[i][j]+D[i-1][j]+D[i][j-1]+D[i-1][j-1]
				a1=0.5*(D[i][j]+D[i][j-1])
				a2=0.5*(D[i-1][j]+D[i][j])
				a3=0.5*(D[i-1][j-1]+D[i-1][j])
				a4=0.5*(D[i][j-1]+D[i-1][j-1])
				R[i][j]=((1/a0)*(a1*V[i+1][j]+a2*V[i][j+1]+a3*V[i-1][j]+a4*V[i][j-1]+rho[i][j]*(h**2)/epsilon_zero)-V[i][j])

def check_if_done():
	done=1
	for i in range(1,l-1):
		for j in range(1,w-1):
			if R[i][j]>error_limit or R[i][j]<-error_limit:
				done=0
	if done==0:
		return False
	else:
		return True

	
def calc_bulk_residuals1D():
	for i in range(l):
		for j in range(1,w-1):
			if BC[i][j]==1:
				R[i][j]=0
				print(i,j)
			else:
				R[i][j]=0.5*(V[i][j+1]+V[i][j-1])-V[i][j]
				
def calc_edge_residuals():
#top edge, to avoid overcounting need to skip the final element in each edge
	for i in range(1,w-1):
		if BC[0][i]==1.:
			R[0][i]=0
		else:
			R[0][i]=0.25*(V[0][i+1]+V[0][i-1]+V[1][i])-V[0][i]
#bottom edge
	for i in range(1,w-1):
		if BC[l-1][i]==1.:
			R[l-1][i]=0
		else:
			R[l-1][i]=0.25*(V[l-1][i+1]+V[l-1][i-1]+V[l-2][i])-V[l-1][i]
#right edge
	for i in range(1,l-1):
		if BC[i][w-1]==1.:
			R[i][w-1]=0
		else:
			R[i][w-1]=0.25*(V[i+1][w-1]+V[i-1][w-1]+V[i][w-2])-V[i][w-1]
#left edge
	for i in range(1,l-1):
		if BC[i][0]==1.:
			R[i][0]=0
		else:
			R[i][0]=0.25*(V[i+1][0]+V[i-1][0]+V[1][1])-V[0][i]
			
#iterate
def iterate(iterations,voltages):
	#must define all voltage plates and charges first
	generate_dirichlet_boundaries()
	for i in range(1,iterations):
		if i/iterations==0.25:
			print("25% complete")
		if i/iterations==.5:
			print("50% complete")
		if i/iterations==.75:
			print("75% complete")
		if i==iterations:
			print("100% complete")
		calc_bulk_residuals()
		#calc_edge_residuals()
		voltages+=R
		
		
def iterate_until_done(voltages):
	#must define all voltage plates and charges first
	generate_dirichlet_boundaries()
	count=20
	for i in range(20):
		calc_bulk_residuals()
		#calc_edge_residuals()
		voltages+=R
	while check_if_done()==False:
		for i in range(1,20):
			calc_bulk_residuals()
			#calc_edge_residuals()
			voltages+=R
			count+=1
	print("This run required ", count," counts to reach the desired accuracy")



#dielectric mesh
D=np.zeros((l,w))
#initialize the dielectric mesh as free space
for i in range(l):
	for j in range(w):
		D[i][j]=1

def create_dielectric_block(x_pos_start, x_length, y_pos_start, y_length,rel_permmitivity):
#this function createsa uniform block of a single permittivity
	for i in range(y_length):
		for j in range(x_length):
			D[i+y_pos_start][j+x_pos_start]=rel_permmitivity

create_dielectric_block(w//4,3*w//4,l//2-5,10,-4)



iterate(100,V)
#iterate_until_done(V)

#electric field calculator

Ex=np.zeros((l,w))
Ey=np.zeros((l,w))
Emagnitude=np.zeros((l,w))

def calc_electric_field():
	for i in range(l-1):
		for j in range(w-1):
			Ex[i][j]=-(V[i][j+1]-V[i][j])/h
			Ey[i][j]=-(V[i+1][j]-V[i][j])/h
			Emagnitude[i][j]=math.sqrt(Ex[i][j]**2+Ey[i][j]**2)


calc_electric_field()

fig, ax = plt.subplots()
# Using matshow here just because it sets the ticks up nicely. imshow is faster.
ax.matshow(Emagnitude, cmap='ocean')
Q=plt.quiver(Ex,Ey)
plt.show()

end=time.clock()
deltatime=end-start
print("this took ", deltatime," s to complete")

import matplotlib.pyplot as plt
import numpy as np


l=10
w=10
grid=np.zeros((l,w))
mass=np.zeros((l,w))

for i in range(l):
	grid[i][(i+1)%l]=1
	grid[i][i]=1
	grid[i][(i-1)%l]=1
	
for i in range(3):
	grid[i+6][i+2]=1
	grid[i+5][i+3]=1
	grid[i+7][i+4]=1
	
print(grid)

#shape=
#	[[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
#[ 0.  0.  0.  0.  1.  0.  0.  0.  0.  0.]
# [ 0.  0.  0.  0.  1.  1.  0.  0.  0.  0.]
 #[ 0.  0.  0.  1.  1.  0.  0.  0.  0.  0.]
 #[ 0.  0.  1.  1.  1.  0.  0.  0.  0.  0.]
# [ 0.  1.  1.  1.  0.  0.  0.  0.  0.  0.]
# [ 1.  1.  1.  1.  0.  0.  0.  0.  0.  0.]
 #[ 1.  1.  1.  1.  1.  1.  1.  1.  1.  1.]
 #[ 1.  1.  1.  1.  1.  1.  1.  1.  1.  1.]

def initialize_grid():
	for i in range(l):
		for j in range(w):
			if grid[i][j]==0:
				print("this is an etchant cell")
				mass[i][j]=1
			elif grid[i][j]==1:
				print("this is a substrate cell")
				mass[i][j]=100
			
def etch():	
	for i in range(l):
		for j in range(w):
			if grid[i][j]==1:
				if grid[(1+i)%l][(j)]==0:
					mass[i][j]-=mass[(i+1)%l][(j)]
				if grid[(i-1)%l][(j)]==0:
					mass[i][j]-=mass[(i-1)%l][(j)]
				if grid[(i)][(j+1)%w]==0:
					mass[i][j]-=mass[(i)][(j+1)%w]
				if grid[(i)][(j-1)%w]==0:
					mass[i][j]-=mass[(i)][(j-1)%w]

def flood():
	for i in range(l):
		for j in range(w):
			if grid[i][j]==1:
				if mass[i][j]<0:
					grid[i][j]=0
					mass[i][j]=1
			
def check_below(a,b):
	if grid[a][b]==1:
		if grid[(a+1)%l][b]==0 and a+1<l:
			mass[a][b]-=mass[a+1][b]
		
def check_above(a,b):
	if grid[a][b]==1:
		if grid[(a-1)%l][b]==0 and a-1>=0:
			mass[a][b]-=mass[a-1][b]
		
def check_right(a,b):
	if grid[a][b]==1:
		if grid[a][(b+1)%w]==0 and b+1<w:
			mass[a][b]-=mass[a][b+1]
		
def check_left(a,b):
	if grid[a][b]==1:
		if grid[a][(b-1)%w]==0 and b-1>=0:
			mass[a][b]-=mass[a][b-1]
			
def etch_edges():
#top edge, to avoid overcounting need to skip the final element in each edge
	for i in range(1,l-1):
		check_right(0,i)
		check_left(0,i)
		check_below(0,i)
		
#bottom edge
	for i in range(1,l-1):
		check_right(l-1,i)
		check_left(l-1,i)
		check_above(l-1,i)
	
#right edge
	for i in range(w):
		check_left(i,w-1)
		check_above(i,w-1)
		check_below(i,w-1)
	
#left edge	
	for i in range(w):
		check_right(i,0)
		check_above(i,0)
		check_below(i,0)
		
def etch_bulk():
	for i in range(l-2):
		for j in range(w-2):
			check_right(i+1,j+1)
			check_left(i+1,j+1)
			check_above(i+1,j+1)
			check_below(i+1,j+1)

initialize_grid()
	
for i in range(50):
	etch_edges()
	etch_bulk()
	flood()

				
print(grid)
print(mass)

fig, ax = plt.subplots()
# Using matshow here just because it sets the ticks up nicely. imshow is faster.
ax.matshow(mass, cmap='seismic')

for (i, j), z in np.ndenumerate(mass):
    ax.text(j, i, '{:0.1f}'.format(z), ha='center', va='center',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='0.7'))

plt.show()

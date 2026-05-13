from mpl_toolkits import mplot3d
import numpy as np 
from math import sqrt
from matplotlib import pyplot as plt 
from JONES import *
f=open('borr_lists.txt','w')


def endpoints(C):
    return C[-1],C[0]
# 0<= t <=1
def intermediate(C,t):
    r1=np.array(endpoints(C)[0])
    r2=np.array(endpoints(C)[1])
    diff=r2-r1
    r=r1+t*diff
    return list(r)

def esets(C):
	ex=np.array([pt[0] for pt in C])
	ey=np.array([pt[1] for pt in C])
	ez=np.array([pt[2] for pt in C])
	return ex, ey, ez

def plotC(C,color):
	ex=esets(C)[0]
	ey=esets(C)[1]
	ez=esets(C)[2]
	ax.plot3D(ex,ey,ez,color)

Col=['r','b','k']

lis=[0,0.1,0.2,0.3,0.4,0.5]
K=np.linspace(0,1,10)
print(K)
LL=[]
for t in K:
    print('t',t)
    # Initial Open Segments
    u1=[[0,0,0],[1,1,0],[2,2,0.5],[3,3,0.5],[4,4,0],[5,5,0],[6,6,0.5],[7,7,0.5],[8,7,0.5],[9,5,0.2],[9,3,0.2],[8,0,0.2],[8,-1,0.2],[6,-1.5,0],[4,-2,0],[2,-1.5,0]]
    u2=[[1,0,0.5],[4,0,0],[5,1,0],[5,4,0.5],[4,5,0.5],[3,6,0],[2,7,0],[-1,6,0],[-1,3,0.5]]
    u3=[[6,0,0.5],[7,6,0],[6,7,0],[3,7,0.5],[2,6,0.5],[2,3,0],[3,2,0],[4,1,0.5]]
    if t>0:
        for l in [u1,u2,u3]:
            l.append(intermediate(l,t))
    BM2=np.array([u1,u2,u3],dtype=object)
    '''PLOT COMMANDS'''
    fig=plt.figure()
    ax=plt.axes(projection='3d')
    ax.set_aspect('equal')
    ax.set_xlim((-4,12))
    ax.set_ylim((-5,15))
    ax.view_init(elev=82,azim=-90)
    ax.set_xlabel("X co-ordinates")
    ax.set_ylabel("Y co-ordinates")
    for i in range(3):
        color=Col[i]
        C1=BM2[i]
        plotC(C1,color)
    fig.savefig('borr_{}.png'.format(round(t,2)))  
    plt.close(fig)
    #JONES POLYNOMIAL CALCULATION
    Object=BM2
    n= 1000

    input_knot=[]

    for arr in Object:
        input_knot.append(list(arr))
    start_time=time.time()
    if t!=K[-1]:
        expected = expected_jones(input_knot, n)

        end_time=time.time()

        '''
        User friendly printing of the polynomial as powers of A vs coeffs
        '''
        LL.append(expected)
print(LL,file=f)
print('',file=f)
print(len(LL),file=f)

f.close()    



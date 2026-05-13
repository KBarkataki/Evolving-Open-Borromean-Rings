import math
import numpy as np
import sys
import random
from copy import copy
import time


# Picks a random point from the surface of a unit sphere.
def get_random_proj():
	b = 0
	while np.sum(b*b)>1 or np.sum(b*b)<1e-4:
		b = 2*np.random.uniform(size=3)-1
	return b/np.sqrt(np.sum(b*b))



# Chooses n = samples points from the surface of a unit sphere.
def fibonacci_sphere(samples):
	points = []
	phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
	for i in range(samples):
		y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
		radius = math.sqrt(1 - y * y)  # radius at y
		theta = phi * i  # golden angle increment
		x = math.cos(theta) * radius
		z = math.sin(theta) * radius
		points.append((x, y, z))
	return points



# Given a normal vector, finds the corresponding 
# plane selecting 2 orthonormal basis vectors representative  of the plane
def get_two_vec(proj_vec): 
	a = np.zeros([3])
	proj_vec=np.array(proj_vec)
	while np.sum(a*a) < 1e-4:
		a = np.random.normal(size=[3])
		a = a - proj_vec*np.sum(a*proj_vec)/np.sum(proj_vec*proj_vec)
	a /= np.sqrt(np.sum(a*a))
	b = np.array([
				  a[1]*proj_vec[2]-a[2]*proj_vec[1], 
				  -a[0]*proj_vec[2]+a[2]*proj_vec[0],
				  a[0]*proj_vec[1]-a[1]*proj_vec[0]
				  ])	#cross product with a and proj
						#orthogonality of a and b can be checked using np.dot(a,b)
	return a, b



#Count how many loops described by the indices
def how_many_loops(inds):
	remaining = np.ones(inds.shape, dtype=np.bool)
	#print('inds',inds)
	#print(remaining)
	count=0
	L=[]
	HL=[]
	while(np.any(remaining)>0):
		next_ind = np.argmax(remaining)
		#print('argmax',np.argmax(remaining))
		length=0
		while(remaining[next_ind]):
			remaining[next_ind] = False
			next_ind = inds[next_ind]
			length+=1
		count+=1
		cyc=[next_ind,next_ind+length]
		L.append(cyc)
	#print('cyc', L)	
	return count, L   

##For open chains
#def open_state_loops(inds.initial):
	#remaining = np.ones(inds.shape, dtype=np.bool)
	#print('inds',inds)
	##print(remaining)
	#count=0
	#L=[]
	#HL=[]
	#while(np.any(remaining)>0):
		#next_ind = np.argmax(remaining)
		#print('argmax',np.argmax(remaining))
		#length=0
		#while(remaining[next_ind]):
			#remaining[next_ind] = False
			#next_ind = inds[next_ind]
			#length+=1
		#count+=1
		#cyc=[next_ind,next_ind+length]
		#L.append(cyc)
	#print('cyc', L)	
	#return count, L       
    



# Reverse the direction of every edge
# i.e. if vertex a went to vertex b, now vertex b goes to vertex a
def reversed_inds(inds):
	return np.arange(inds.shape[0])[np.argsort(inds)]



#This generates a matrix indicating which edges are crossing, where they are crossing, 
#which line is on top, etc.
''' CRAMER'S RULE 2x2'''
def det2x2(A):
	assert A.shape == (2,2)
	return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def Cramer(A):
	assert A.shape == (2,3)
	D = det2x2(A[:,:2])
	if D == 0:
		return
	Dx = det2x2(A[:,[2,1]])
	Dy = det2x2(A[:,[0,2]])
	return Dx*1.0/D, Dy*1.0/D



# Reidemeister MOVES

def RM1(bool_mask,inds):
	if np.any(bool_mask):
		edge1 = np.argmax(np.any(bool_mask, 0))
		edge2 = np.argmax(bool_mask[edge1,:])
		Row1=bool_mask[inds[edge1] : , inds[edge2] :]
		#print('Row1',Row1)
		dummy=1
		#for index in range(len(Row1)):
		if np.any(Row1):
#			if Row1[index]==True and index!= edge2:
#				dummy=0
			dummy=0
			#print('dummy',dummy)
		if dummy!=0:
			#print('yes')
			bool_mask[edge1,edge2]=False
			bool_mask[edge2,edge1]=False
	#print('bool_mask',bool_mask)		
	return bool_mask

def RM2(bool_mask, over_or_under,inds):
	if np.any(bool_mask):
		edge1 = np.argmax(np.any(bool_mask, 0))
		edge2 = np.argmax(bool_mask[edge1,:])
		Temp=bool_mask[inds[edge1] : , inds[edge2] :]
		if np.any(Temp):
			tm3=np.argmax(np.any(Temp, 0))
			edge3 = inds[edge1]+ np.argmax(np.any(Temp, 0))
			edge4 = inds[edge2]+np.argmax(Temp[tm3,:])
			if over_or_under[edge1,edge2]==over_or_under[edge3,edge4]:
				Row1=bool_mask[inds[edge1] : edge3 , inds[edge2] : edge4]
				dummy=1
				if np.any(Row1):
					dummy=0
				if dummy!=0:
					bool_mask[edge1,edge2]=False
					bool_mask[edge2,edge1]=False
					bool_mask[edge3,edge4]=False
					bool_mask[edge4,edge3]=False
		#print('bool_mask',bool_mask)		
	return bool_mask

def simplification(BM,over_or_under,inds):
	for k in range(2):
		for i in range(10):
			BM=RM1(BM,inds)
#		for j in range(10):
#			BM=RM2(BM,over_or_under,inds)
	return BM		


def max_len(Ch):
	l=0
	for ch in Ch:
		l=max([l,len(ch)])
	return l

def J_mult(P,Q):
	Z={}
	for i in P:
		for j in Q: 
			try:
				Z[i+j]+=P[i]*Q[j]
			except:
				Z[i+j]=P[i]*Q[j]
	return Z

def J_add(P,Q):
	Z=P
	for j in Q: 
		try:
			Z[j]+=Q[j]
		except:
			Z[j]=Q[j]
	return Z

def J_smult(a,P):
	Z={}
	for i in P:
		Z[i]=a*P[i]
	return Z
	
def dfactor(N):
	dpoly={0:1}
	for i in range(N):
		dpoly=J_mult(dpoly,{-2:-1,2:-1})
	return dpoly



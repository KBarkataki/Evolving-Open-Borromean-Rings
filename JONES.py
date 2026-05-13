'''
Given an entangled object with finite componenets in 3-space the following program
calculates its JONES POLYNOMIAL
'''
import math
import numpy as np
import sys
import random
from copy import copy
import time
#from collections import Counter
from functions import *


closed=0 # if open, then assign value 0
''' 
================================================================
JONES POLYNOMIAL CALCULATION

In expected_Jones(Ch,n),
Ch (short for chain) signifies a 
list (of different components) of lists (of co-ordinates in any particular component) 
n signifies the no. of projections being considered (for averaging out our result)

In get_jones_poly(coords,proj_vector,inds),
coords signifies the co-ordinates of points on the curve in 3space
proj_vector signifies a particular direction of projection
inds signifies the connections among the points in the curve in concern.
'''

def expected_jones(Ch, n):

	#initial connections among INDICES corresponding to the co-ords in 3-space
	inds=np.array([])
	ini=0
	p=np.concatenate(Ch)
	for ch in Ch:
		ch=np.array(ch)
		sub_inds=np.concatenate([np.arange(1+ini, ch.shape[0]+ini), np.array([ini])])
		inds=np.concatenate([inds,sub_inds])
		ini=len(inds)	 
	inds=np.array([int(i) for i in inds])	

	# Uniformly GENERATE n vectors on the unit sphere
	#if n>1:
		#points=fibonacci_sphere(n)
##		points=[]
##		for i in range(n):
##			points.append(get_random_proj())
	if n>1:
		#points=fibonacci_sphere(n)
		points=[]
		for i in range(n):
			theta=random.uniform(0,2*math.pi)
			'''range of phi is adjusted math.pi/6 was original choice'''			
			phi=random.uniform(0,math.pi/2)
			xcoord=math.cos(phi)*math.cos(theta)
			ycoord=math.cos(phi)*math.sin(theta)
			zcoord=math.sin(phi)
			points.append(np.array([xcoord,ycoord,zcoord]))
	else: # choose one at random
		points=[np.array([0.0,0.0,1.0])]
		#points=[get_random_proj()]
	#print('points',points)
#	base = np.array([0])


	JPOLY=[]
	Jnone=0
	for proj_vector in points:
		#print('PRO')
		JPoly=[]
		rj = get_jones_poly(p, proj_vector, inds) # JP along a particular proj_vector
		if rj!=None:
			for base in rj:
				lenn=int((len(base)-1)/2)
				powers=range(-lenn,lenn+1,1)
				poly={}
				for i in range(len(base)):
					if base[i]!=0:
						poly[powers[i]]=base[i]
				JPoly.append(poly)
			JPoly.append(dfactor(len(JPoly)-1))
			#print('Jpoly',JPoly)
			Zpoly={0:1}
			for i in range(len(JPoly)):
				Zpoly=J_mult(JPoly[i],Zpoly)
			JPOLY.append(Zpoly)
		else:
			Jnone+=1

	Zpoly={0:0}
	for i in range(len(JPOLY)):
		Zpoly=J_add(JPOLY[i],Zpoly)
	n=n-Jnone
	#print('Jnone',Jnone,n)
	Zpoly=J_smult(1./n,Zpoly)
	return Zpoly
	
		#print("Jones Polynomial : ", Zpoly)

#		if base.shape[0] < rj.shape[0]:
#			pad_width = int((rj.shape[0] - base.shape[0])/2)
#			base = np.concatenate([np.zeros(pad_width), base, np.zeros(pad_width)])
#		elif base.shape[0] > rj.shape[0]:
#			pad_width = int((base.shape[0] - rj.shape[0])/2)
#			rj = np.concatenate([np.zeros(pad_width), rj, np.zeros(pad_width)])
#		base += rj
	#return base/n


def get_jones_poly(coords, proj_vector, inds):
	# plane for a particular KNOT DIAGRAM
	x, y = get_two_vec(proj_vector) 
	#proj gives the projection of the 3D coords into the x-y plane
	proj = np.matmul(coords, np.stack([x, y], 1))
	# CHECK if the projection is valid in the sense of KNOT DIAGRAM.
	Check_proj=np.unique([str(i) for i in proj],return_counts=True)
	print('Check_proj',Check_proj)
	#print(max(Check_proj[1]))
	print(max(Check_proj[1]))
	if max(Check_proj[1])>1:
		return None
	# We store depth (along new z vector i.e. proj_vector) of each vertex 
	depth = np.matmul(coords, np.expand_dims(proj_vector, 1))[:,0]	
	bool_mask, over_or_under, u, right_or_left, proj, depth, inds = get_bool_overlap_etc(proj, inds, depth,start=0)

	#print("before", bool_mask)
	#bool_mask=simplification(bool_mask,over_or_under,inds)
	#print(np.count_nonzero(bool_mask)/2)
	if  np.count_nonzero(bool_mask)>40:
		return None
	else:
	#''' DISCONNECTED COMPONENTS IN LINKS/LINKOIDS'''
		partition=how_many_loops(inds)[1]
		#print("PART",partition)
		conn_dict={}
		for s in range(len(partition)):
			conn_dict[s]=[s]
			BMs=bool_mask[partition[s][0] : partition[s][1], :]
			for t in range(len(partition)):
				if t!=s:
					col_i=partition[t][0]
					col_f=partition[t][1]
					if np.any(BMs[:,col_i : col_f]):
						conn_dict[s].append(t)
		#print('DICT',conn_dict)
		disco=[]
		discard=[]
		for i in conn_dict:
			if i not in discard:
				list_i=[]
				for j in range(i,len(conn_dict)):
					if len(set(conn_dict[i]).intersection(set(conn_dict[j]))) > 0:
						list_i=set(list_i).union(set(conn_dict[j]))
				discard=discard+list(list_i)
				disco.append(list(list_i))
		#print('DISCO',disco)
		#print('DISCO',disco, max_len(disco))	
		
		if max_len(disco)<=20:	# scope for tweaking
			K_list=[]
			for compo in disco[:1]:
				#print('compo',compo)
				compo_proj=[]
				compo_depth=[]
				for chn in compo:
					tup=partition[chn]
					compo_proj.append(list(proj[tup[0]:tup[1]]))
					#print('DP',len(depth),len(depth[tup[0]:tup[1]]))
					compo_depth+=list(depth[tup[0]:tup[1]])
				compo_inds=np.array([])
				ini=0
				#print('compo_proj',compo_proj)
				p=np.concatenate(compo_proj)
				for chn in compo_proj:
					chn=np.array(chn)
					compo_sub_inds=np.concatenate([np.arange(1+ini, chn.shape[0]+ini), np.array([ini])])
					compo_inds=np.concatenate([compo_inds,compo_sub_inds])
					ini=len(compo_inds)	 
				compo_inds=np.array([int(i) for i in compo_inds])	
				#print('lens',len(compo_inds),len(p),len(compo_depth))
				#print('inds',compo_inds)
			
				bool_mask, over_or_under, u, right_or_left, Proj, Depth, inds = get_bool_overlap_etc(proj, inds, depth,start=0)	
			
				a = get_partial_poly(bool_mask, over_or_under, right_or_left, inds)
				#print('a',a)
				# Writhe calculation
				b = get_writhe(bool_mask,over_or_under,right_or_left)
				print('b',b)
				if b > 0:
					without_quarter = np.concatenate([np.zeros(6*b), a*(-1)**b])
					K= without_quarter
				elif b < 0:
					without_quarter = np.concatenate([a*(-1)**(-b), np.zeros(-6*b)])
					K= without_quarter
				else:
					K= a
				K_list.append(K)
			#print("KKK",K_list)
			return(K_list)
		else:
			return None
		#return A

''' WRITHE '''
def get_writhe(bool_mask,over_or_under,right_or_left):
	bool_mask[np.arange(bool_mask.shape[0]).reshape((-1,1)) <= np.arange(bool_mask.shape[0]).reshape((1, -1))] = False
	return np.sum(2*np.int32(over_or_under == right_or_left)[bool_mask]-1)

'''BRACKET POLYNOMIAL'''
'''
Gets the characteristic polynomial
proj: should be a nx2 matrix of vertices
inds: describes the graph through the indices (should be a closed knot) (n size array)
depth: the depth of each vertex in three space (important for crossings)
edge2ignore: this value should indicate what edge is dropped (preventing a cycle)
'''
def get_partial_poly(bool_mask, over_or_under, right_or_left, inds):

	'''
	Check if there are any intersections
	'''
	#print(bool_mask.shape)
	if np.any(bool_mask):
		#print('IN YES')
		edge1 = np.argmax(np.any(bool_mask, 0))
		edge2 = np.argmax(bool_mask[edge1,:])
		#print("edge1: ", edge1, " edge2: ", edge2)

		bool_mask1 = np.copy(bool_mask)
		bool_mask1[edge1, edge2] = False
		bool_mask1[edge2, edge1] = False

		'''
		We are going to create two new paths: path inds1 and path inds2
		'''
		inds1 = np.copy(inds)
		inds2 = np.copy(inds)
		reversed_ind = reversed_inds(inds)
		#print("reversed_ind", reversed_ind)

		#inds1 is easy, because we just swap the destinations of edge1 and edge2
		inds1[edge1] = inds[edge2]
		inds1[edge2] = inds[edge1]

		#inds2 is substantially trickier
		inds2[edge1] = edge2
		first2flip = edge2
		#We have to reverse a few of the directions (hence a while loop)
		replacement = np.arange(inds.shape[0])
		#print('replacement',replacement)
		replacement[edge1] = edge2
		while first2flip != inds[edge1] and first2flip != inds[edge2]:
			inds2[first2flip] = reversed_ind[first2flip]
			replacement[first2flip] = reversed_ind[first2flip]

			first2flip = reversed_ind[first2flip]
		
		if first2flip == inds[edge1]:
			inds2[inds[edge1]] = inds[edge2]
			replacement[inds[edge1]] = edge1
		else:
			inds2[inds[edge2]] = inds[edge1]
			replacement[inds[edge2]] = edge1
			
		bool_mask2 = bool_mask1[replacement[:,None], replacement[None,:]]
		#print('bool 2',bool_mask2)
		right_or_left2 = right_or_left[replacement[:,None], replacement[None,:]]
		need_to_swap = replacement!=np.arange(inds.shape[0])
		swap_mask = need_to_swap[:,None] != need_to_swap[None,:]
		right_or_left2[swap_mask] = np.logical_not(right_or_left2[swap_mask])
		over_or_under2 = over_or_under[replacement[:,None], replacement[None,:]]
		#We need the sub-polynomials
		partial_poly1 = get_partial_poly(bool_mask1, np.copy(over_or_under), np.copy(right_or_left), inds1)
		partial_poly2 = get_partial_poly(bool_mask2, over_or_under2, right_or_left2, inds2)
		#print('PP1',partial_poly1,'PP2',partial_poly2)
		#This part adds the polynomials
		#We first need to make sure all the polynomials are of the same degree
		#We do this by padding with zeros when necessary
		if partial_poly1.shape[0] > partial_poly2.shape[0]:
			width = int((partial_poly1.shape[0] - partial_poly2.shape[0])/2)
			partial_poly2 = np.concatenate([np.zeros([width]), partial_poly2, np.zeros([width])])
		if partial_poly2.shape[0] > partial_poly1.shape[0]:
			width = int((partial_poly2.shape[0] - partial_poly1.shape[0])/2)
			partial_poly1 = np.concatenate([np.zeros([width]), partial_poly1, np.zeros([width])])

		#We check for rightedness with aboveness

		if over_or_under[edge1, edge2] == right_or_left[edge1, edge2]:
			return (
				np.concatenate([partial_poly1, np.zeros([2])]) + 
				np.concatenate([np.zeros([2]), partial_poly2])
			)
		else:
			return (
				np.concatenate([partial_poly2, np.zeros([2])]) + 
				np.concatenate([np.zeros([2]), partial_poly1])
			)
	else:
		base = np.array([1])
		#We just count the number of loops. We subtract by 1 because one of the 
		#loops is broken up by the missing link
		num_loops = how_many_loops(inds)[0] - 1 # FOR CLOSED LOOPS
		#print(num_loops)
		for powi in range(num_loops):
			base = np.concatenate([-base, np.zeros([4])]) + np.concatenate([np.zeros([4]), -base])
		return base


def get_bool_overlap_etc(proj, inds, depth,start=0):
	#print('depth along proj axis : ', depth )
	# directed initial to final 
	#print('IN')
	#print('lens',len(proj),len(depth),len(inds))
	dif = proj[inds] - proj 
	'''
	bool_mask is a nxn matrix: the value at row i, column j indicates 
	if edge i intersects with edge j
	'''
	bool_mask = np.zeros((len(inds),len(inds)), dtype=bool)
	'''
	If there is an intersection, which edge is above the other?
	'''
	over_or_under = np.zeros((len(inds),len(inds)), dtype=bool)
	'''
	u matrix stores intersection points between edges
	'''
	u = np.zeros((len(inds),len(inds)), dtype=object)
	#print("dif[:,1] ",dif[:,1])
	#print("dif[:,0] ",dif[:,0])
	flip_neg = np.stack([dif[:,1], -dif[:,0]], 1)
	#print('flip_neg', flip_neg) 
	# Finding intersection point b/w edge i and edge j
	for i in range(len(inds)):
		for j in range(inds[i],len(inds),1):
			pi=proj[i][0]*proj[inds[i]][1]-proj[i][1]*proj[inds[i]][0]
			pj=proj[j][0]*proj[inds[j]][1]-proj[j][1]*proj[inds[j]][0]
			Aij=np.array([[flip_neg[i][0],flip_neg[i][1],pi],[flip_neg[j][0],flip_neg[j][1],pj]])
			result=Cramer(Aij)
			#print(i,j,result)
			if result != None:
				# parametrisation of intersection point
				# either proj[inds[i]][0]-proj[i][0] is nonzero or proj[inds[i]][1]-proj[i][1] is nonzero
				xdif_i=proj[inds[i]][0]-proj[i][0]
				xdif_j=proj[inds[j]][0]-proj[j][0]
				if xdif_i!=0:
					t=(result[0]-proj[i][0])/(proj[inds[i]][0]-proj[i][0])
				else:
					t=(result[1]-proj[i][1])/(proj[inds[i]][1]-proj[i][1])
				if xdif_j!=0:
					s=(result[0]-proj[j][0])/(proj[inds[j]][0]-proj[j][0])
				else:
					s=(result[1]-proj[j][1])/(proj[inds[j]][1]-proj[j][1])
				zero_val=1.e-6 # buffer required for numerical process
				if zero_val<t<1-zero_val and zero_val<s<1-zero_val: 
					bool_mask[i,j]=True
					bool_mask[j,i]=True
					zi=depth[i]+t*(depth[inds[i]]-depth[i])
					zj=depth[j]+s*(depth[inds[j]]-depth[j])
					if zi>zj:
						over_or_under[i,j]=True # edge i over edge j
					else:
						over_or_under[j,i]=True # edge j over edge i
					u[i,j]=result
					u[j,i]=result
	'''
	We need to know which edge is to the right of the other
	'''
	right_or_left = np.sum(np.expand_dims(flip_neg, 0)*np.expand_dims(dif, 1), -1)>0

	# This loop runs only for Open curves. Edits crossings 
	if closed==0:	
		partition=how_many_loops(inds)[1]
		edges_2_ignore=[tup[1]-1 for tup in partition] # false edges
		
		for i in edges_2_ignore:
			bool_mask[i,:]=False
			bool_mask[:,i]=False
		#print('edges_2_ig',edges_2_ignore)
	else:
		partition=how_many_loops(inds)[1]
		edges_2_ignore=[]

	#print('PART',partition)

	''' FIX POLY '''
	# Introduction of new vertices wherever there is 1 edge crossing more than 1 edge.
	if start ==0:
		dim=bool_mask.shape[0]

		new_vectors = []
		new_depth = []
		#Check if there are any edge with >=2 intersections
		for row in range(dim):
			new_vectors.append(list(proj[row]))
			new_depth.append(depth[row])
			if row not in edges_2_ignore:
				pts_on_segment=[]	# before ordering
				for col in range(dim):
					if col not in edges_2_ignore:
						if bool_mask[row,col]==True:
							pts_on_segment.append(u[row,col])
				

				if len(pts_on_segment)>1:
					Chk_proj=np.unique([str(i) for i in proj],return_counts=True)
					if max(Chk_proj[1])>1:
						return None
					#print('edge, int',row, pts_on_segment)
					#print('YES')
					xdif_0=proj[inds[row]][0]-proj[row][0]
					xdif_1=proj[inds[row]][1]-proj[row][1]
					mag_vs_pt={}	# magnitude of displacement vs point
					for elt in pts_on_segment:
						mag_vs_pt[abs(elt[0]-proj[row][0])]=elt
					mag_keys=list(mag_vs_pt.keys())
					mag_keys.sort()
					Ordered_points=[]
					for mag in mag_keys:
						Ordered_points.append(mag_vs_pt[mag])
					for i in range(len(Ordered_points)-1):
						p1=Ordered_points[i]
						p2=Ordered_points[i+1]
						new_p=[0,0]
						new_p[0]=(p1[0]+p2[0])/2.
						new_p[1]=(p1[1]+p2[1])/2.
						if xdif_0!=0:
							new_dep=depth[row]+(new_p[0]-proj[row][0])/xdif_0*(depth[inds[row]]-depth[row])
						else:
							new_dep=depth[row]+(new_p[1]-proj[row][1])/xdif_1*(depth[inds[row]]-depth[row])

						new_vectors.append(new_p)
						new_depth.append(new_dep)
		#print('number of elts in the set new_vectors',len(new_vectors))
		NNV=[]
		#print('part',partition)
		for chunk in partition:
			#print('chunk',chunk)
			first_pt=list(proj[chunk[0]])
			#print('first pount',first_pt)
			pt_ind=new_vectors.index(first_pt)
			#print('index of first point',pt_ind)
			NNV.append(pt_ind)
			#print('each step NV',NV)
		NNV.append(len(new_vectors))

		''' CAUTION. Apparently there are repititions in NV. Hence this segment of code'''
		NV=[]
		for x in NNV:
			if x not in NV:
				NV.append(x)

		#print('all NVS',NV)
		inds=np.array([])
		for i in range(len(NV)-1):
			sub_inds=np.concatenate([np.arange(NV[i]+1, NV[i+1]), np.array([NV[i]])])
			#print('SUB',sub_inds)
			inds=np.concatenate([inds,sub_inds])
			#print('----',inds)
		#print('before array', len(inds))	 
		inds=np.array([int(i) for i in inds])
		Check_freq=np.unique(inds,return_counts=True)
		#print('CFF',Check_freq,max(Check_freq[0]),max(Check_freq[1]))
		#print('index511', list(inds).index(511))
		#print('after array', inds, len(inds))	
		proj=np.array(new_vectors)
		depth=np.array(new_depth)
		#print('nnn',len(inds),len(depth),len(proj))

		bool_mask, over_or_under, u, right_or_left, proj, depth, inds = get_bool_overlap_etc(proj, inds, depth,start=1)

	return bool_mask, over_or_under, u, right_or_left, proj, depth, inds



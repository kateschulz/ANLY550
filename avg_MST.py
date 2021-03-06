#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 16:14:39 2019

@author: kateschulz
"""

import sys
from math import sqrt
import itertools
import random
from functools import lru_cache


#Compute 2D euclidean distance
def two_d_dist(x1,y1,x2,y2):
     dist = sqrt((x1 - x2)**2 + (y1-y2)**2) 
     return dist

#Compute 3D euclidean distance
def three_d_dist(x1,y1,z1,x2,y2,z2):
     dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
     return dist
 
#Compute 4D euclidean distance
def four_d_dist(x1,y1,z1,t1,x2,y2,z2, t2):
     dist = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2 + (t2 - t1)**2)
     return dist  

#Compute prime factors of number n
def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors 

# Create adjacency matrix for prims algorithm taking a set of vertices V and
# graph G with edge weights
def create_adj_matrix(V, G):
  
  adj_matrix = []
  
  # create n x n matrix filled with 0 edge weights between vertices
  for i in range(0, V):
    adj_matrix.append([])
    for j in range(0, V):
      adj_matrix[i].append(0)
      
  # populate matrix with computed edge weights
  for i in range(0, len(G)):
    adj_matrix[G[i][0]][G[i][1]] = G[i][2]
    adj_matrix[G[i][1]][G[i][0]] = G[i][2]
    
  return adj_matrix

# Create graph with list of edge weights as input
def create_graph(edges):
    #create n distinct vertices for edge weights
    new_vertices = random.sample(range(0, n), n) 

    #create tuples of vertices for graph
    new_vertex_pairs = list(itertools.combinations(new_vertices,2))
    list1, list2 = zip(*new_vertex_pairs)
    
    # create graph from vertices and edge weights
    tuples = list(zip(list1,list2,edges))
    graph = [list(i) for i in tuples]
    
    return graph

############## CACHE create_edge_weights ###############
@lru_cache(maxsize=10000000)

# compute edge weigths based on n and dim
def create_edge_weights(n, dim):
    # Case dim = 0
    if dim == 0:
        #create n distinct vertices 
        vertices = random.sample(range(0, n), n) 
        
        #create pairs of vertices for complete graph
        vertex_pairs = list(itertools.combinations(vertices,2))
        list1, list2 = zip(*vertex_pairs)
        
        # n choose 2 edge weights from [0,1]
        edge_weights = []
        for j in range(len(vertex_pairs)):
            edge_weights.append(random.uniform(0, 1))
                
        
    elif (dim == 2 or dim == 3 or dim == 4):
        #create n*dim points from [0,1]
        vertex_points = []
        for j in range(n*dim):
            vertex_points.append(random.uniform(0, 1))

        # pair up points and combine pairs for dim
        vertex_iter = [iter(vertex_points)] * dim 
        vertices = zip(*vertex_iter)
        vertex_pairs = list(itertools.combinations(vertices, 2))
            
        # Case dim = 2
        if dim == 2:
            #compute 2D euclidean distance
            edge_weights = []
            for j in range(len(vertex_pairs)):
                edge_weights.append(two_d_dist(vertex_pairs[j][0][0],vertex_pairs[j][0][1],
                                               vertex_pairs[j][1][0],vertex_pairs[j][1][1]))
            
        # Case dim = 3
        elif dim == 3: 
            #compute 3D euclidean distance
            edge_weights = []
            for j in range(len(vertex_pairs)):
                edge_weights.append(three_d_dist(vertex_pairs[j][0][0],vertex_pairs[j][0][1],
                                                 vertex_pairs[j][0][2], vertex_pairs[j][1][0],
                                                 vertex_pairs[j][1][1], vertex_pairs[j][1][2]))
        # Case dim = 4
        elif dim == 4:
            #compute 4D euclidean distance
            edge_weights = []
            for j in range(len(vertex_pairs)):
                edge_weights.append(four_d_dist(vertex_pairs[j][0][0],vertex_pairs[j][0][1],
                                                vertex_pairs[j][0][2], vertex_pairs[j][0][3],
                                                vertex_pairs[j][1][0],vertex_pairs[j][1][1], 
                                                vertex_pairs[j][1][2], vertex_pairs[j][1][3]))
                
        # Only accept dim input of 0, 2, 3, or 4
        else:
            return str("Please enter a dimension of 0, 2, 3, or 4!")
    
    return edge_weights
 
# Run Prim's algorithm on a set of V vertices and graph G with vertex pairs
# and edge weights
def prims(V, G):
  
  # create adj matrix from graph
  adj_matrix = create_adj_matrix(V, G)
  
  # choose initial vertex from graph
  vertex = 0
  
  # initialize MST, edges, and track visited vertices and min edge
  MST = []
  edges = []
  visited = []
  min_edge = [None,None,float('inf')]
  
  # run prims algorithm until we create an MST
  # that contains every vertex from the graph
  while len(MST) != V-1:
    
    # mark vertex as visited
    visited.append(vertex)
    
    # add each edge to list of potential edges in MST
    for r in range(0, V):
      if adj_matrix[vertex][r] != 0:
        edges.append([vertex,r,adj_matrix[vertex][r]])
        
    # find edge with the smallest weight to a vertex that is not visited 
    for e in range(0, len(edges)):
      if edges[e][2] < min_edge[2] and edges[e][1] not in visited:
        min_edge = edges[e]
        
    # remove min weight edge from list of edges
    edges.remove(min_edge)

    # add min edge to MST
    MST.append(min_edge)
      
    # start at new vertex and reset min edge
    vertex = min_edge[1]
    min_edge = [None,None,float('inf')]
  
  return MST
      
############## CACHE avg_MST ###############
@lru_cache(maxsize=10000000)

# Run Prim's algorithm to compute the MST weight from n vertices in a given 
# dimension over a given number of trials
def avg_MST(n, dim, num_trials):
    
    #compute average weight over trials
    weight = []
    for i in range(num_trials):
    #compute edge weights
        edge_weights = create_edge_weights(n, dim)
        
        # for n sufficiently large, prune the number of edges
        if n >= (2 ** (7)):
            twos = prime_factors(n)
            m = len(twos)
            k_n = 0.5 ** (m-6)
            
            for item in edge_weights:
                if item >= k_n:
                    edge_weights.remove(item)
    
        # create graph         
        graph = create_graph(edge_weights)

        # run Prim's algorithm and compute MST weight
        path_weight = prims(n, graph)
        for i in range(len(path_weight)):
            weight.append(path_weight[i][2])
        
    #Compute average weight over all trials
    return sum(weight)/num_trials
    
# Read std.in parameters
n = int(sys.argv[2])
num_trials = int(sys.argv[3])
dim = int(sys.argv[4])

#print average MST weight and given parameters
print(avg_MST(n,dim,num_trials), n, num_trials, dim)



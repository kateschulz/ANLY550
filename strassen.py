#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:09:50 2019

@author: kateschulz
"""

import sys 
from math import ceil, floor, log, sqrt

#function to read input .txt file
def read(filename, dim):
    lines = open(filename, 'r').read().splitlines()
    
    #first n^2 lines are A, second n^2 lines are B
    A = lines[0:int(dim)**2]
    A = [float(i) for i in A]
    
    B = lines[int(dim)**2:2*int(dim)**2]
    B = [float(i) for i in B]
   
    return A, B


#function to take dimension n and flat list and create nxn matrix
def matrix_create(dim, data_list):
    mat = []
    for i in range(dim):
        rows = []
        for j in range(dim):
            rows.append(data_list[dim * i + j])
        mat.append(rows)

    return mat
 
#function to multiply two nxn matrices      
def matrix_mult(A, B): 
    n = len(A)
    C = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for k in range(n):
            for j in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

#function to add two nxn matrices
def matrix_add(A, B):
    n = len(A)
    C = [[0 for j in range(0, n)] for i in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            C[i][j] = A[i][j] + B[i][j]
    return C

#function to subtract two nxn matrices
def matrix_subtract(A, B):
    n = len(A)
    C = [[0 for j in range(0, n)] for i in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            C[i][j] = A[i][j] - B[i][j]
    return C
 
#optimized Strassen's algorithm with cross-over point n_0      
def strassen_opt(A, B, n_opt):
    
    n_0 = n_opt 
    n = len(A)

    # below cross-over point, perform naive matrix multiplication
    if n <= n_0:
        return matrix_mult(A, B)
    
    # above cross-over point, perform Strassen's algorithm
    else:
        # divide A and B into submatrices
        sub_dim = n//2
        a11 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        a12 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        a21 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        a22 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]

        b11 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        b12 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        b21 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        b22 = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]

        a_sub = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        b_sub = [[0 for j in range(0, sub_dim)] for i in range(0, sub_dim)]
        
        for i in range(0, sub_dim):
            for j in range(0, sub_dim):
                a11[i][j] = A[i][j]            
                a12[i][j] = A[i][j + sub_dim]    
                a21[i][j] = A[i + sub_dim][j]  
                a22[i][j] = A[i + sub_dim][j + sub_dim] 

                b11[i][j] = B[i][j]           
                b12[i][j] = B[i][j + sub_dim]   
                b21[i][j] = B[i + sub_dim][j]    
                b22[i][j] = B[i + sub_dim][j + sub_dim] 
            
        #compute p1 through p7 using Strassen's equations:
        a_sub = matrix_add(a11, a22)
        b_sub = matrix_add(b11, b22)
        p1 = strassen_opt(a_sub, b_sub)

        a_sub = matrix_add(a21, a22)
        p2 = strassen_opt(a_sub, b11)

        b_sub = matrix_subtract(b12, b22) 
        p3 = strassen_opt(a11, b_sub) 

        b_sub = matrix_subtract(b21, b11) 
        p4 =strassen_opt(a22, b_sub)  

        a_sub = matrix_add(a11, a12)     
        p5 = strassen_opt(a_sub, b22)     

        a_sub = matrix_subtract(a21, a11) 
        b_sub = matrix_add(b11, b12)    
        p6 = strassen_opt(a_sub, b_sub) 

        a_sub = matrix_subtract(a12, a22) 
        b_sub = matrix_add(b21, b22)     
        p7 = strassen_opt(a_sub, b_sub) 
        
        # compute product C of sub_matrices of A and B
        c12 = matrix_add(p3, p5) 
        c21 = matrix_add(p2, p4) 

        a_sub = matrix_add(p1, p4)
        b_sub = matrix_add(a_sub, p7)
        c11 = matrix_subtract(b_sub, p5)

        a_sub = matrix_add(p1, p3) 
        b_sub = matrix_add(a_sub, p6)
        c22 = matrix_subtract(b_sub, p2) 

        C = [[0 for j in range(0, n)] for i in range(0, n)]
        for i in range(0, sub_dim):
            for j in range(0, sub_dim):
                C[i][j] = c11[i][j]
                C[i][j + sub_dim] = c12[i][j]
                C[i + sub_dim][j] = c21[i][j]
                C[i + sub_dim][j + sub_dim] = c22[i][j]
        return C

#function to pad matrices before calling optimized Strassen's
def strassen(A, B):
    n = len(A) 
    
    #handle powers of 2 cases with no padding
    if ((ceil(log(n,2)) == floor(log(n,2)))):  
        # empirical n_0 is 32 
        C = strassen_opt(A,B, 32)
        return C
    
    #handle all other cases by padding with zeros
    else:
        #pad to increase dimensions to next power of 2
        power_two = lambda n: 2**int(ceil(log(n,2)))
        m = power_two(n)

        A_pad = [[0 for i in range(m)] for j in range(m)]
        B_pad = [[0 for i in range(m)] for j in range(m)]
        
        for i in range(n):
            for j in range(n):
                A_pad[i][j] = A[i][j]
                B_pad[i][j] = B[i][j]
        # empirical n_0 is 64
        C_pad = strassen_opt(A_pad, B_pad, 64)

    #extract C from C padded with zeros
    C = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = C_pad[i][j]
    return C
    
#function to print diagonal of matrix
def print_diagonal(matrix):
    for i in range(len(matrix)):
        print (matrix[i][i])
    
#run program in shell with user input        
a, b = read(sys.argv[3], sys.argv[2])    
  
#set up matrices 
A = matrix_create(int(sqrt(len(a))),a)
B = matrix_create(int(sqrt(len(b))),b)

#run optimized strassen's algorithm
C = strassen(A, B)    
print_diagonal(C)







                

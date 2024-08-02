#Ioannis Kavoukis, 6318
#Geometry for Beam Crane

import numpy as np
import AEM6318_AEM_data_calculator as dataFun

def massEq(AEM):
    [_,_,_,_,_,_,_,A_bar,_,n,_] = dataFun.AEM_dataCalc(AEM)
    
    a = (10 + 10*np.sqrt(2)/3 + 2*np.sqrt(3)/3 + 2*np.sqrt(5) + 2*np.sqrt(6) + 
         8*n + 10*np.sqrt(2)/3*n)/(8 + 2*np.sqrt(5) + 2*np.sqrt(6) + 2*np.sqrt(2) + 8*n)
    
    A_beam = a*A_bar    
    
    return A_beam

def craneGeometry(AEM):

    [F,E,A,B,phi,theta,A_diag,A_bar,L,n,v] = dataFun.AEM_dataCalc(AEM)
    A_beam = massEq(AEM)    
    
    #NODES
    base_nodes = 2
    support_nodes = 3
    top_nodes = 3
    mainBody_nodes = (n+1)*4
    n_nodes = base_nodes + support_nodes + top_nodes + mainBody_nodes + 1
    
    nodes_id = np.array(range(1,n_nodes+1))
    nodes_coord = np.zeros((n_nodes,6)) #coordinates
    nodes_BC = np.zeros((n_nodes,6)) #boundary conditions
    nodes_F = np.zeros((n_nodes,6)) #forces Fx Fy Fz
    nodes = list((nodes_id,nodes_coord,nodes_BC,nodes_F))

    #Base
    nodes_coord[0,:] = 0 #the center of my coordinate system
    nodes_coord[1,:] = 0
    nodes_coord[1,1] = L

    #Main Body
    i = 0
    for i in range(n+1):
        
        nodes_coord[2+4*i:6+4*i,2] = L*(i+1) #z coordinate
        
        nodes_coord[2+i*4,0] = -L/2 #x coordinate
        nodes_coord[3+i*4,0] = -L/2
        nodes_coord[4+i*4,0] = L/2
        nodes_coord[5+i*4,0] = L/2

        nodes_coord[2+i*4,1] = L #y coordinate
        nodes_coord[3+i*4,1] = 0
        nodes_coord[4+i*4,1] = 0
        nodes_coord[5+i*4,1] = L
        
    
    #Top
    nodes_coord[1+(n+1)*4+1:1+(n+1)*4+3,0] = L/2
    nodes_coord[1+(n+1)*4+3,0] = L/2+L
    
    nodes_coord[1+(n+1)*4+1,1] = L
    nodes_coord[1+(n+1)*4+2,1] = 0
    nodes_coord[1+(n+1)*4+3,1] = L/2

    nodes_coord[1+(n+1)*4+1:1+(n+1)*4+3,2] = L*(n+2)
    nodes_coord[1+(n+1)*4+3,2] = L*(n+2)-L/2

    #Support
    nodes_coord[n_nodes-4:n_nodes-1,0] = -A
    
    nodes_coord[n_nodes-4,1] = 0
    nodes_coord[n_nodes-3,1] = L/2
    nodes_coord[n_nodes-2,1] = L
    
    nodes_coord[n_nodes-4:n_nodes-1,2] = B
    
    #Coordination point    
    nodes_coord[n_nodes-1,0:2] = L/4
    nodes_coord[n_nodes-1,2] = 5*L
    
    #Rotation
    j = 0
    for j in range(2,n_nodes-4):
        X1 = nodes_coord[j,0]*np.cos(np.pi/2-phi) + nodes_coord[j,2]*np.sin(np.pi/2-phi)
        Z1 = nodes_coord[j,2]*np.cos(np.pi/2-phi) - nodes_coord[j,0]*np.sin(np.pi/2-phi)
        nodes_coord[j,0] = X1
        nodes_coord[j,2] = Z1
            
    k = 0
    for k in range(1,n_nodes-1):
        X2 = nodes_coord[k,0]*np.cos(theta) + nodes_coord[k,1]*np.sin(-theta)
        Y2 = nodes_coord[k,1]*np.cos(theta) - nodes_coord[k,0]*np.sin(-theta)
        nodes_coord[k,0] = X2
        nodes_coord[k,1] = Y2
    
    
    #ELEMENTS
    base_elements = 5
    support_elements = 4
    top_elements = 13
    mainBody_elements = n*8
    n_elements = base_elements + support_elements + top_elements + mainBody_elements
    
    elements_id = np.array(range(1,n_elements+1))
    elements_nodes = np.zeros((n_elements,3))        
    elements_area = np.zeros((n_elements,1))
    elements_elasticity = E*np.ones((n_elements,1))
    elements_poisson = v*np.ones((n_elements,1))

    elements = list((elements_id , elements_nodes , elements_area , elements_elasticity , elements_poisson))
    
    #Base
    elements_nodes[0,0] = 1
    elements_nodes[0,1] = 2
    
    elements_nodes[1,0] = 1
    elements_nodes[1,1] = 4
    
    elements_nodes[2,0] = 1
    elements_nodes[2,1] = 5
    
    elements_nodes[3,0] = 2
    elements_nodes[3,1] = 6
    
    elements_nodes[4,0] = 2
    elements_nodes[4,1] = 3
    
    #Main Body
    i = 0; j = 0
    
    for i in range(n):
        k = 0; m=0
        for j in range(3,7):
            elements_nodes[base_elements+k+8*i,0] = j+4*i
            elements_nodes[base_elements+k+8*i,1] = (j+1)+4*i
            if j == 6:
                elements_nodes[base_elements+k+8*i,1] = (j-3)+4*i
                
            elements_nodes[base_elements+4+k+8*i,0] = j+4*i
            elements_nodes[base_elements+4+k+8*i,1] = (j+4)+4*i

            k = k + 1

    #Top
    i = 0; j = 0; k = 0; m = 0
    
    for i in range(3,8):
        elements_nodes[base_elements+mainBody_elements+k,0] = i+4*n
        elements_nodes[base_elements+mainBody_elements+k,1] = i+4*n+1
        
        k = k + 1

    for j in range(3,6,2):
        elements_nodes[base_elements+mainBody_elements+5+m,0] = j+4*n
        elements_nodes[base_elements+mainBody_elements+5+m,1] = j+4*n+3
        
        m = m + 1
        
    j = 0; m = 0    
    for j in range(5,9):
        elements_nodes[base_elements+mainBody_elements+7+m,0] = j+4*n
        elements_nodes[base_elements+mainBody_elements+7+m,1] = n_nodes-4
                
        m = m + 1
        
    j = 0; m = 0    
    for j in range(3,5):
        elements_nodes[base_elements+mainBody_elements+11+m,0] = j+4*n
        elements_nodes[base_elements+mainBody_elements+11+m,1] = j+4*n+4
         
        m = m + 1

    #Support
    lvl = np.floor(0.4*n) #level on which the 2 'orange' supports will adhere
    
    elements_nodes[n_elements-4,0] = 3+4*lvl
    elements_nodes[n_elements-4,1] = n_nodes-1
    
    elements_nodes[n_elements-3,0] = 4+4*lvl
    elements_nodes[n_elements-3,1] = n_nodes-3
    
    elements_nodes[n_elements-2,0] = 3+4*n
    elements_nodes[n_elements-1,0] = 4+4*n  
    elements_nodes[n_elements-2:n_elements,1] = n_nodes-2
    
    
    elements_area[0:n_elements-4] = A_beam
    elements_area[n_elements-4:n_elements] = A_bar

    i = 0
    for i in range(n_elements):
        elements_nodes[i,2] = n_nodes

    #Default Forces
    nodes_F[n_nodes-5,2] = -F
    
    #Default Boundary Conditions
    nodes[2][0,0:3] = 1
    nodes[2][1,0:3] = 1
    nodes[2][1,1] = 0
    nodes[2][n_nodes-4:n_nodes-1,0:3] = 1
    
    
    return nodes , elements
    
    
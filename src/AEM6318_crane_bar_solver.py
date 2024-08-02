#Ioannis Kavoukis, 6318
#Bar Crane Solver

import numpy as np

# Import data from txt files
nodes_file = open('text_files/AEM6318_crane_bar_Nodes.txt','r')

nFile = nodes_file.read().splitlines()
n_nodes = len(nFile)

nodes_id = np.array(range(1,n_nodes+1))
nodes_coord = np.zeros((n_nodes,3))
nodes_BC = np.zeros((n_nodes,3))
nodes_F = np.zeros((n_nodes,3))

i = 0
for i in range(n_nodes):
    dummyVar = nFile[i].split(",")
    # nodes_id[i] = int(dummyVar[0])
    
    nodes_cord = np.fromstring(dummyVar[1][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(3):
        nodes_coord[i][j] = nodes_cord[j]
        
    nodes_bc = np.fromstring(dummyVar[2][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(3):
        nodes_BC[i][j] = nodes_bc[j]
           
    nodes_forces = np.fromstring(dummyVar[3][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(3):
        nodes_F[i][j] = nodes_forces[j]

nodes_file.close
nodes = list((nodes_id,nodes_coord,nodes_BC,nodes_F))


elements_file = open('text_files/AEM6318_crane_bar_Elements.txt','r')

eFile = elements_file.read().splitlines()
n_elements = len(eFile)

elements_id = np.array(range(1,n_elements+1))
elements_nodes = np.zeros((n_elements,2))        
elements_area = np.zeros((n_elements,1))
elements_elasticity = np.zeros((n_elements,1))

i = 0
for i in range(n_elements):
    dummyVar = eFile[i].split(",")
    # elements_id[i] = int(dummyVar[0])
    
    el_nodes = np.fromstring(dummyVar[1][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(2):
        elements_nodes[i][j] = el_nodes[j]
        
    elements_area[i] = float(dummyVar[2][1:-1])
    elements_elasticity[i] = float(dummyVar[3][1:-1])

elements_file.close
elements = list((elements_id , elements_nodes , elements_area , elements_elasticity))

###############################################################################

#Solver
K_stiff = np.zeros((len(nodes_id)*3,len(nodes_id)*3))

elements_length = np.zeros((len(elements_id)))
i = 0
for i in range(len(elements_id)):
    elements_length[i] = np.linalg.norm(nodes_coord[int(elements_nodes[i,1])-1,:]-nodes_coord[int(elements_nodes[i,0])-1,:])
    
    x_length = (nodes_coord[int(elements_nodes[i,0])-1,0]-nodes_coord[int(elements_nodes[i,1])-1,0])/elements_length[i]
    y_length = (nodes_coord[int(elements_nodes[i,0])-1,1]-nodes_coord[int(elements_nodes[i,1])-1,1])/elements_length[i]
    z_length = (nodes_coord[int(elements_nodes[i,0])-1,2]-nodes_coord[int(elements_nodes[i,1])-1,2])/elements_length[i]

    ke = (elements_area[i]*elements_elasticity[i]/elements_length[i])*np.array([[1 , -1],[-1 , 1]])
    T = np.array([[x_length,y_length,z_length,0,0,0],[0,0,0,x_length,y_length,z_length]])
    Ke = np.transpose(T)@ke@T

    j = 3*(int(elements_nodes[i,0])-1)
    k = 3*(int(elements_nodes[i,1])-1)
    
    K_stiff[j:j+3,j:j+3] = K_stiff[j:j+3,j:j+3] + Ke[0:3,0:3]
    K_stiff[j:j+3,k:k+3] = K_stiff[j:j+3,k:k+3] + Ke[0:3,3:6]
    K_stiff[k:k+3,j:j+3] = K_stiff[k:k+3,j:j+3] + Ke[3:6,0:3]
    K_stiff[k:k+3,k:k+3] = K_stiff[k:k+3,k:k+3] + Ke[3:6,3:6]

#Penalty Method, where there is BC, multiply K(i,j) with large number
K = K_stiff.copy()
i = 0; j = 0
for i in range(len(nodes_id)):
    for j in range(3):
        if nodes_BC[i,j] == 1:
            K[3*i+j,3*i+j] += 1e7 #check if indices are correct

forces = np.reshape(nodes_F, [nodes_F.size,1])

D = np.linalg.det(K)
if D == 0:
    print("The K matrix is not invertible")

u = np.linalg.inv(K)@forces
U = np.reshape(u, [-1,3])

Rs = K_stiff@u

strain = np.zeros(len(elements_id))
stress = np.zeros(len(elements_id))

nodes_disp_coord = nodes_coord + U

elements_def_length = np.zeros(len(elements_id))
i = 0
for i in range(len(elements_id)):
    elements_def_length[i] = np.linalg.norm(nodes_disp_coord[int(elements_nodes[i,1]-1),:]-nodes_disp_coord[int(elements_nodes[i,0]-1),:])
    
    strain[i] = (elements_def_length[i]-elements_length[i])/elements_length[i]
    stress[i] = elements_elasticity[i]*strain[i]

reactions = np.reshape(Rs,[-1,3])

reactions_total = np.zeros(len(reactions))
i = 0
for i in range(len(reactions)):
    reactions_total[i] = np.linalg.norm(reactions[i,:])

###############################################################################

#output of solver
strains_file = open('text_files/AEM6318_crane_bar_strains.txt','w+')
i = 0
for i in range(len(strain)):
    strains_file.write(""+str(strain[i])+"\n")
strains_file.close()

stress_file = open('text_files/AEM6318_crane_bar_stress.txt','w+')
i = 0
for i in range(len(stress)):
    stress_file.write(""+str(stress[i])+"\n")
stress_file.close()

displacement_file = open('text_files/AEM6318_crane_bar_displacement.txt','w+')
i = 0
for i in range(len(U)):
    displacement_file.write(""+str(U[i,:])+"\n")
displacement_file.close()

reactions_file = open('text_files/AEM6318_crane_bar_reactions.txt','w+')
i = 0
for i in range(len(reactions)):
    reactions_file.write(""+str(reactions[i,0])+" "+str(reactions[i,1])+" "+str(reactions[i,2])+" "+str(reactions_total[i])+"\n")
reactions_file.close()





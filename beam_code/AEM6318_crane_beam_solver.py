#Ioannis Kavoukis, 6318
#Beam Crane Solver

import numpy as np

# Import data from txt files
nodes_file = open('text_files/AEM6318_crane_beam_Nodes.txt','r')
nFile = nodes_file.read().splitlines()
n_nodes = len(nFile)

nodes_id = np.array(range(1,n_nodes+1))
nodes_coord = np.zeros((n_nodes,6))
nodes_BC = np.zeros((n_nodes,6))
nodes_F = np.zeros((n_nodes,6))

i = 0
for i in range(n_nodes):
    dummyVar = nFile[i].split(",")
    # nodes_id[i] = int(dummyVar[0])
    
    nodes_cord = np.fromstring(dummyVar[1][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(6):
        nodes_coord[i][j] = nodes_cord[j]
        
    nodes_bc = np.fromstring(dummyVar[2][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(6):
        nodes_BC[i][j] = nodes_bc[j]
           
    nodes_forces = np.fromstring(dummyVar[3][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(6):
        nodes_F[i][j] = nodes_forces[j]

nodes_file.close
nodes = list((nodes_id,nodes_coord,nodes_BC,nodes_F))


elements_file = open('text_files/AEM6318_crane_beam_Elements.txt','r')
eFile = elements_file.read().splitlines()
n_elements = len(eFile)

elements_id = np.array(range(1,n_elements+1))
elements_nodes = np.zeros((n_elements,3))        
elements_area = np.zeros((n_elements,1))
elements_elasticity = np.zeros((n_elements,1))
elements_poisson = np.zeros((n_elements,1))

i = 0
for i in range(n_elements):
    dummyVar = eFile[i].split(",")
    # elements_id[i] = int(dummyVar[0])
    
    el_nodes = np.fromstring(dummyVar[1][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(3):
        elements_nodes[i][j] = el_nodes[j]
        
    elements_area[i] = float(dummyVar[2][1:-1])
    elements_elasticity[i] = float(dummyVar[3][1:-1])
    elements_poisson[i] = float(dummyVar[4][1:-1])

elements_file.close
elements = list((elements_id , elements_nodes , elements_area , elements_elasticity , elements_poisson))

###############################################################################
 
#Solver
K_stiff = np.zeros(((len(nodes[0])-1)*6,(len(nodes[0])-1)*6))

elements_length = np.zeros((len(elements[0])))
i = 0
for i in range(len(elements[0])-4):
    elements_length[i] = np.linalg.norm(nodes[1][int(elements[1][i,1])-1,:]-nodes[1][int(elements[1][i,0])-1,:])
   
    X1 = nodes[1][int(elements[1][i,0])-1,0]; X2 = nodes[1][int(elements[1][i,1])-1,0]; X3 = nodes[1][int(elements[1][i,2])-1,0]
    Y1 = nodes[1][int(elements[1][i,0])-1,1]; Y2 = nodes[1][int(elements[1][i,1])-1,1]; Y3 = nodes[1][int(elements[1][i,2])-1,1]
    Z1 = nodes[1][int(elements[1][i,0])-1,2]; Z2 = nodes[1][int(elements[1][i,1])-1,2]; Z3 = nodes[1][int(elements[1][i,2])-1,2]

    A123 = np.sqrt(((Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1))**2+((Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1))**2+((X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1))**2)

    lx = (X2-X1)/elements_length[i]
    mx = (Y2-Y1)/elements_length[i]
    nx = (Z2-Z1)/elements_length[i]
    
    lz = ((Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1))/A123
    mz = ((Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1))/A123
    nz = ((X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1))/A123
    
    ly = mz*nx-nz*mx
    my = nz*lx-lz*nx
    ny = lz*mx-mz*lx

    T3 = np.array([[lx , mx , nx],[ly , my , ny],[lz , mz , nz]])
    T = np.block([[T3, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
                  [np.zeros((3,3)), T3, np.zeros((3,3)), np.zeros((3,3))],
                  [np.zeros((3,3)), np.zeros((3,3)), T3, np.zeros((3,3))],
                  [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), T3]])

    #basic characteristics
    a = elements_length[i]/2
    r = np.sqrt(elements[2][i]/np.pi)
    Iy = np.pi*r**4/4 #moment of inertia at y axis
    Iz = Iy  #same moment of inertia because of cylindrical cross section
    J = 2*Iy #polar moment of inertia
    E = elements[3][i]
    v = elements[4][i]
    G = E/(2*(1+v))
    
    #helpful constants
    c1 = (elements[2][i]*E)/(2*a)
    c2 = (3*E*Iz)/(2*a**3)
    c3 = (3*E*Iz)/(2*a**2)
    c4 = (3*E*Iy)/(2*a**3)
    c5 = (3*E*Iy)/(2*a**2)
    c6 = (G*J)/(2*a)
    c7 = (E*Iy)/a
    c8 = (E*Iz)/a

    ke = np.array([[c1, 0, 0, 0, 0, 0, -c1, 0, 0, 0, 0, 0],
                   [0, c2, 0, 0, 0, c3, 0, -c2, 0, 0, 0, c3],
                   [0, 0, c4, 0, -c5, 0, 0, 0, -c4, 0, -c5, 0],
                   [0, 0, 0, c6, 0, 0, 0, 0, 0, -c6, 0, 0],
                   [0, 0, -c5, 0, 2*c7, 0, 0, 0, c5, 0, c7, 0],
                   [0, c3, 0, 0, 0, 2*c8, 0, -c3, 0, 0, 0, c8],
                   [-c1, 0, 0, 0, 0, 0, c1, 0, 0, 0, 0, 0],
                   [0, -c2, 0, 0, 0, -c3, 0, c2, 0, 0, 0, -c3],
                   [0, 0, -c4, 0, c5, 0, 0, 0, c4, 0, c5, 0],
                   [0, 0, 0, -c6, 0, 0, 0, 0, 0, c6, 0, 0],
                   [0, 0, -c5, 0, c7, 0, 0, 0, c5, 0, 2*c7, 0],
                   [0, c3, 0, 0, 0, c8, 0, -c3, 0, 0, 0, 2*c8]],dtype=object)

    Ke = np.transpose(T)@ke@T

    j = 6*(int(elements[1][i,0])-1)
    k = 6*(int(elements[1][i,1])-1)
    K_stiff[j:j+6,j:j+6] = K_stiff[j:j+6,j:j+6] + Ke[0:6,0:6]
    K_stiff[j:j+6,k:k+6] = K_stiff[j:j+6,k:k+6] + Ke[0:6,6:12]
    K_stiff[k:k+6,j:j+6] = K_stiff[k:k+6,j:j+6] + Ke[6:12,0:6]
    K_stiff[k:k+6,k:k+6] = K_stiff[k:k+6,k:k+6] + Ke[6:12,6:12]

i = 0
for i in range(len(elements[0])-4,len(elements[0])):
    elements_length[i] = np.linalg.norm(nodes[1][int(elements[1][i,1])-1,:]-nodes[1][int(elements[1][i,0])-1,:])
   
    X1 = nodes[1][int(elements[1][i,0])-1,0]; X2 = nodes[1][int(elements[1][i,1])-1,0]; #X3 = nodes[1][int(elements[1][i,2])-1,0]
    Y1 = nodes[1][int(elements[1][i,0])-1,1]; Y2 = nodes[1][int(elements[1][i,1])-1,1]; #Y3 = nodes[1][int(elements[1][i,2])-1,1]
    Z1 = nodes[1][int(elements[1][i,0])-1,2]; Z2 = nodes[1][int(elements[1][i,1])-1,2]; #Z3 = nodes[1][int(elements[1][i,2])-1,2]

    l = (X2-X1)/elements_length[i]
    m = (Y2-Y1)/elements_length[i] 
    n = (Z2-Z1)/elements_length[i] 

    kk = elements[3][i]*elements[2][i]/elements_length[i]
    mat1 = np.array([[l**2, l*m, l*n],[l*m, m**2, m*n],[l*n, m*n, n**2]])
    O3 = np.zeros((3,3))
    I3 = np.eye(3)
    
    Ke = kk*np.block([[mat1,O3,-mat1,O3],[O3,I3,O3,O3],[-mat1,O3,mat1,O3],[O3,O3,O3,I3]])

    j = 6*(int(elements[1][i,0])-1)
    k = 6*(int(elements[1][i,1])-1)
    K_stiff[j:j+6,j:j+6] = K_stiff[j:j+6,j:j+6] + Ke[0:6,0:6]
    K_stiff[j:j+6,k:k+6] = K_stiff[j:j+6,k:k+6] + Ke[0:6,6:12]
    K_stiff[k:k+6,j:j+6] = K_stiff[k:k+6,j:j+6] + Ke[6:12,0:6]
    K_stiff[k:k+6,k:k+6] = K_stiff[k:k+6,k:k+6] + Ke[6:12,6:12]

#Penalty Method, where there is BC, multiply K(i,j) with large number
K = K_stiff.copy()
i = 0; j = 0
for i in range(len(nodes[0])-1):
    for j in range(6):
        if nodes[2][i,j] == 1:
            K[6*i+j,6*i+j] += 1e7 

forces = np.reshape(nodes[3][0:len(nodes[3])-1,:], [nodes[3].size-6,1])

D = np.linalg.det(K)
if D == 0:
    print("The K matrix is not invertible")

u = np.linalg.inv(K)@forces

U = np.reshape(u, [-1,6])
U = np.block([[U],[np.zeros((1,6))]])
Rs = K_stiff@u

strain = np.zeros(len(elements[0]))
stress = np.zeros(len(elements[0]))

nodes_disp_coord = nodes[1] + U

elements_def_length = np.zeros(len(elements[0]))
i = 0
for i in range(len(elements[0])):
    elements_def_length[i] = np.linalg.norm(nodes_disp_coord[int(elements[1][i,1]-1),:]-nodes_disp_coord[int(elements[1][i,0]-1),:])
    
    strain[i] = (elements_def_length[i]-elements_length[i])/elements_length[i]
    stress[i] = elements[3][i]*strain[i]

reactions = np.reshape(Rs,[-1,6])

reactions_f_total = np.zeros(len(reactions))
reactions_m_total = np.zeros(len(reactions))
i = 0
for i in range(len(reactions)):
    reactions_f_total[i] = np.linalg.norm(reactions[i,0:3])
    reactions_m_total[i] = np.linalg.norm(reactions[i,3:6])

###############################################################################

#output of solver
strains_file = open('text_files/AEM6318_crane_beam_strains.txt','w+')
i = 0
for i in range(len(strain)):
    strains_file.write(""+str(strain[i])+"\n")
strains_file.close()

stress_file = open('text_files/AEM6318_crane_beam_stress.txt','w+')
i = 0
for i in range(len(stress)):
    stress_file.write(""+str(stress[i])+"\n")
stress_file.close()

displacement_file = open('text_files/AEM6318_crane_beam_displacement.txt','w+')
i = 0
for i in range(len(U)):
    disp = str(U[i,:])
    disp = disp.replace("\n","")
    displacement_file.write(""+disp+"\n")
displacement_file.close()

reactions_file = open('text_files/AEM6318_crane_beam_reactions.txt','w+')
i = 0
for i in range(len(reactions)):
    reactions_file.write(""+str(reactions[i,0])+" "+str(reactions[i,1])+" "+str(reactions[i,2])+" \
                          "+str(reactions[i,3])+" "+str(reactions[i,4])+" "+str(reactions[i,5])+" \
                          "+str(reactions_f_total[i])+" "+str(reactions_m_total[i])+"\n")
reactions_file.write("0  0  0  0  0  0  0  0\n")
reactions_file.close()




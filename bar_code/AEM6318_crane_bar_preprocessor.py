#Ioannis Kavoukis, 6318
#Bar Crane Pre-Processor

import AEM6318_AEM_data_calculator as dataFun
import AEM6318_crane_bar_geometry as geometryFun

import math

def numel(n):
    if n > 0:
        digits = int(math.log10(n))+1
    elif n == 0:
        digits = 1
    else:
        digits = int(math.log10(-n))+1

    return digits

a = 0;
while a != 1:
    AEM = int(input("Give me your AEM:\n") or "6318")
    if numel(int(AEM)) == 4:
        a = 1

[F,E,A,B,phi,theta,A_diag,A_bar,L,n] = dataFun.AEM_dataCalc(AEM)

#give me the geometry as input
nodes , elements = geometryFun.craneGeometry(AEM)
n_nodes = len(nodes[0])
n_elements = len(elements[0])

###############################################################################

#PLOTS
import matplotlib.pyplot as plt
i = 0
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(len(elements[0])):
    if elements[2][i] == A_bar:
        ax.plot(xs=[nodes[1][int(elements[1][i,0]-1),0], nodes[1][int(elements[1][i,1]-1),0]], \
                ys=[nodes[1][int(elements[1][i,0]-1),1], nodes[1][int(elements[1][i,1]-1),1]], \
                zs=[nodes[1][int(elements[1][i,0]-1),2], nodes[1][int(elements[1][i,1]-1),2]], color = 'blue', linewidth = 1.0)

    elif elements[2][i] == A_diag:
        ax.plot(xs=[nodes[1][int(elements[1][i,0]-1),0], nodes[1][int(elements[1][i,1]-1),0]], \
                ys=[nodes[1][int(elements[1][i,0]-1),1], nodes[1][int(elements[1][i,1]-1),1]], \
                zs=[nodes[1][int(elements[1][i,0]-1),2], nodes[1][int(elements[1][i,1]-1),2]], color = 'red', linewidth = .5)

j = 0       
for j in range(n_nodes):
    ax.scatter(nodes[1][j,0],nodes[1][j,1],nodes[1][j,2],s=5, color = 'black')

plt.xlim(-0, 20000)
plt.ylim(-1000, 20000)    
ax.view_init(25, -45)
plt.xlabel("x")
plt.ylabel("y")
plt.savefig('a.pdf',format="pdf")
plt.show()

###############################################################################
#IF I USE THIS, THEN COMMENT THE DEFAULT SETTINGS IN AEM6318_crane_bar_geometry
# #Forces
# numF = int(input("Give me the number of Applied Forces:\n") or "1")

# i = 0
# for i in range(numF):
#     fVar = str(i+1)
#     a = 0
#     while a != 1:
#         dirF = input("Choose the direction of the Applied Force ("+fVar+")\n(in x|y|z):  ")  or "z"
#         if dirF == "x":
#             a = 1
#             forceDir = 0
#         elif dirF == "y":
#             a = 1
#             forceDir = 1
#         elif dirF == "z":
#             a = 1
#             forceDir = 2
 
#     a = 0
#     while a != 1:
#         forceNode = -1 + int(input("Choose the node of the Applied Force ("+fVar+")\n(in 1-"+str(n_nodes)+"):  ") or " "+str(n_nodes-3)+" ") 
#         if 0 <= forceNode <= n_nodes-1:
#             a = 1
            
#     forceMagn = int(input("Give me the magnitude of the Applied Force ("+fVar+")\n(in Newton):  ") or "-20000")
    
#     nodes[3][forceNode,forceDir] = forceMagn #construct node-force matrix


# #Boundary Conditions
# numBC_nodes = int(input("Give me the number of nodes with Boundary Conditions:\n") or "5")

# i = 0 
# for i in range(numBC_nodes):
#     bcVar = str(i+1)
#     a = 0
#     while a != 1:
#         bcNode = -1 + int(input("Choose the node of the Boundary Condition ("+bcVar+")\n(in 1-"+str(n_nodes)+"):  "))
#         if 0 <= bcNode <= n_nodes-1:
#             a = 1
        
#     numBC = int(input("How many DoF do you want to restrict for node "+str(bcNode+1)+"\n(in 1-3):  "))
#     if numBC == 3:
#         nodes[2][bcNode,:] = 1 #construct node-BC matrix
#     else:
#         j = 0
#         for j in range(numBC):
#             dirBC = input("Choose a direction you want to restrict for node "+str(bcNode+1)+"\n(in x|y|z):  ")
#             if dirBC == "x":
#                 a = 1
#                 bcDir = 0
#             elif dirBC == "y":
#                 a = 1
#                 bcDir = 1
#             elif dirBC == "z":
#                 a = 1
#                 bcDir = 2
    
#             nodes[2][bcNode,bcDir] = 1 #construct node-BC matrix

###############################################################################

# Produce the pre-processor output
import os

if os.path.exists("text_files") == False:
    os.makedirs("text_files")
    
nodes_file = open('text_files/AEM6318_crane_bar_Nodes.txt','w+')
i = 0
for i in range(n_nodes):
    node_id = str(nodes[0][i])
    node_cord = str(nodes[1][i,:])
    node_BC = str(nodes[2][i,:])
    node_force = str(nodes[3][i,:])
    nodes_file.write(""+node_id+","+node_cord+","+node_BC+","+node_force+"\n")
nodes_file.close()
  
elements_file = open('text_files/AEM6318_crane_bar_Elements.txt','w+')
i = 0
for i in range(n_elements):
    elements_id = str(elements[0][i])
    elements_cord = str(elements[1][i,:])
    elements_area = str(elements[2][i])
    elements_elast = str(elements[3][i])
    elements_file.write(""+elements_id+","+elements_cord+","+elements_area+","+elements_elast+"\n")
elements_file.close()









        
        
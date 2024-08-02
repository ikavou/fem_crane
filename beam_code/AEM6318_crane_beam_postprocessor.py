#Ioannis Kavoukis, 6318
#Beam Crane Post-Processor

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Import data from txt files
scale = float(input("Use a scale factor for the graphs:\n") or "1")

nodes_file = open('text_files/AEM6318_crane_beam_Nodes.txt','r')
disp_file = open('text_files/AEM6318_crane_beam_displacement.txt','r')
reactions_file = open('text_files/AEM6318_crane_beam_reactions.txt','r')
nFile = nodes_file.read().splitlines()
dFile = disp_file.read().splitlines()
rFile = reactions_file.read().splitlines()
n_nodes = len(nFile)

nodes_id = np.array(range(1,n_nodes+1))
nodes_coord = np.zeros((n_nodes,6))
nodes_BC = np.zeros((n_nodes,6))
nodes_F = np.zeros((n_nodes,6))
nodes_disp_coord_scaled = np.zeros((n_nodes,6))
U = np.zeros((n_nodes,6))
nodes_reactions = np.zeros((n_nodes,8))

i = 0
for i in range(n_nodes):
    dummyVar = nFile[i].split(",")

    nodes_id[i] = int(dummyVar[0])
    
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
        
    displacement = np.fromstring(dFile[i][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(6):
        U[i][j] = displacement[j]
        nodes_disp_coord_scaled[i][j] = nodes_coord[i][j] + scale*displacement[j]

    react = np.fromstring(rFile[i][0:], dtype=float, sep=' ')
    j = 0
    for j in range(8):
        nodes_reactions[i][j] = react[j]

nodes_file.close
nodes = list((nodes_id,nodes_coord,nodes_disp_coord_scaled,nodes_BC,nodes_F,nodes_reactions))


elements_file = open('text_files/AEM6318_crane_beam_Elements.txt','r')
strain_file = open('text_files/AEM6318_crane_beam_strains.txt','r')
stress_file = open('text_files/AEM6318_crane_beam_stress.txt','r')
eFile = elements_file.read().splitlines()
strainFile = strain_file.read().splitlines()
stressFile = stress_file.read().splitlines()
n_elements = len(eFile)

elements_id = np.array(range(1,n_elements+1))
elements_nodes = np.zeros((n_elements,3))        
elements_area = np.zeros((n_elements,1))
elements_elasticity = np.zeros((n_elements,1))
elements_poisson = np.zeros((n_elements,1))
elements_strain = np.zeros((n_elements,1))
elements_stress = np.zeros((n_elements,1))

i = 0
for i in range(n_elements):
    dummyVar = eFile[i].split(",")

    elements_id[i] = int(dummyVar[0])
    
    el_nodes = np.fromstring(dummyVar[1][1:-1], dtype=float, sep=' ')
    j = 0
    for j in range(3):
        elements_nodes[i][j] = el_nodes[j]
        
    elements_area[i] = float(dummyVar[2][1:len(dummyVar[2])-1])
    elements_elasticity[i] = float(dummyVar[3][1:len(dummyVar[3])-1])
    elements_poisson[i] = float(dummyVar[4][1:len(dummyVar[4])-1])
    elements_strain[i] = float(strainFile[i])
    elements_stress[i] = float(stressFile[i])

elements_file.close
strain_file.close()
stress_file.close()
elements = list((elements_id, elements_nodes, elements_area, elements_elasticity, elements_poisson, elements_strain, elements_stress))

###############################################################################

#PLOTS

###################Undeformed_Deformed_Crane###################################
# i = 0; j = 0
# fig1 = plt.figure()
# ax = fig1.add_subplot(111, projection='3d')
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# for i in range(len(elements[0])):
#     ax.plot(xs=[nodes[1][int(elements[1][i,0]-1),0], nodes[1][int(elements[1][i,1]-1),0]], \
#             ys=[nodes[1][int(elements[1][i,0]-1),1], nodes[1][int(elements[1][i,1]-1),1]], \
#             zs=[nodes[1][int(elements[1][i,0]-1),2], nodes[1][int(elements[1][i,1]-1),2]], color = 'royalblue', linewidth = 0.5)

#     ax.plot(xs=[nodes[2][int(elements[1][i,0]-1),0], nodes[2][int(elements[1][i,1]-1),0]], \
#             ys=[nodes[2][int(elements[1][i,0]-1),1], nodes[2][int(elements[1][i,1]-1),1]], \
#             zs=[nodes[2][int(elements[1][i,0]-1),2], nodes[2][int(elements[1][i,1]-1),2]], color = 'blue', linewidth = 1)
        
# j = 0
# for j in range(n_nodes-1):
#     ax.scatter(nodes[2][j,0],nodes[2][j,1],nodes[2][j,2],s=4, color = 'black')
#     # ax.scatter(nodes[1][j,0],nodes[1][j,1],nodes[1][j,2],s=4, color = 'grey')
    
# ax.view_init(25 , -40)
# plt.xlabel("x")
# plt.ylabel("y")

# plt.title("Undeformed-Deformed Structure\n(scale="+str(scale)+")")
# plt.savefig('Undeformed_Deformed_Beam_Crane.pdf',format="pdf")
# plt.show()

# #########################Forces-Reactions######################################
fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
plt.xlim(-3000, 18000)
plt.ylim(-3000, 18000)
i = 0
for i in range(len(elements[0])):
    ax.plot(xs=[nodes[2][int(elements[1][i,0]-1),0], nodes[2][int(elements[1][i,1]-1),0]], \
            ys=[nodes[2][int(elements[1][i,0]-1),1], nodes[2][int(elements[1][i,1]-1),1]], \
            zs=[nodes[2][int(elements[1][i,0]-1),2], nodes[2][int(elements[1][i,1]-1),2]], color = 'blue', linewidth = 1)
 
j = 0
for j in range(n_nodes-1):
    ax.scatter(nodes[2][j,0],nodes[2][j,1],nodes[2][j,2],s=4, color = 'black')

i = 0
for i in range(n_nodes-1):
    if (nodes[4][i]).any != 0:
        coordF = nodes[2][i,:]
        ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.3, color = 'black')

    if abs(nodes[5][i,6]) > 1e-5:
        coordR = nodes[2][i,:]
        ax.quiver(coordR[0],coordR[1],coordR[2],nodes[5][i,0],nodes[5][i,1],nodes[5][i,2], length=5e-7*nodes[5][i,6], color = 'red')

#CHECK SCALE OF ARROWS
ax.view_init(25 , -50)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Forces & Reactions")
plt.savefig('Forces_Reactions_Beam_Crane.pdf',format="pdf")
plt.show()

# #############################Strains###########################################
# fig3 = plt.figure()
# ax = fig3.add_subplot(111, projection='3d')

# i = 0
# for i in range(len(elements[0])):
#     xs = (nodes[2][int(elements[1][i,0]-1),0], nodes[2][int(elements[1][i,1]-1),0])
#     ys = (nodes[2][int(elements[1][i,0]-1),1], nodes[2][int(elements[1][i,1]-1),1])
#     zs = (nodes[2][int(elements[1][i,0]-1),2], nodes[2][int(elements[1][i,1]-1),2])

#     strainRange = (elements[5][i]-min(elements[5]))/(max(elements[5])-min(elements[5]))
#     ax.plot(xs, ys, zs, color = plt.cm.jet(strainRange), linewidth = 1)

# i = 0
# for i in range(n_nodes-1):
#     if (nodes[4][i]).any != 0:
#         coordF = nodes[2][i,:]
#         ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.2, color = 'black')

# j = 0
# for j in range(n_nodes-1):
#     ax.scatter(nodes[2][j,0],nodes[2][j,1],nodes[2][j,2],s=4, color = 'black')

# norm = mpl.colors.Normalize(min(elements[5]), max(elements[5]))
# cmap = plt.get_cmap('jet', 100)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig3.colorbar(sm, ax=ax, location = 'left')
# cbar.set_label('\u03B5 [mm/mm]', rotation=90)
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# plt.title("Axial Strain")
# ax.view_init(25 , -50)
# plt.xlabel("x")
# plt.ylabel("y")


# plt.savefig('Strains_Beam_Crane.pdf',format="pdf")
# plt.show()

# ################################Stress#########################################
# fig4 = plt.figure()
# ax = fig4.add_subplot(111, projection='3d')

# i = 0
# for i in range(len(elements[0])):
#     xs = (nodes[2][int(elements[1][i,0]-1),0], nodes[2][int(elements[1][i,1]-1),0])
#     ys = (nodes[2][int(elements[1][i,0]-1),1], nodes[2][int(elements[1][i,1]-1),1])
#     zs = (nodes[2][int(elements[1][i,0]-1),2], nodes[2][int(elements[1][i,1]-1),2])

#     stressRange = (elements[6][i]-min(elements[6]))/(max(elements[6])-min(elements[6]))
#     ax.plot(xs, ys, zs, color = plt.cm.jet(stressRange), linewidth = 1)

# j = 0
# for j in range(n_nodes-1):
#     ax.scatter(nodes[2][j,0],nodes[2][j,1],nodes[2][j,2],s=4, color = 'black')

# i = 0
# for i in range(n_nodes-1):
#     if (nodes[4][i]).any != 0:
#         coordF = nodes[2][i,:]
#         ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.2, color = 'black')

# norm = mpl.colors.Normalize(min(elements[6]), max(elements[6]))
# cmap = plt.get_cmap('jet', 100)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig4.colorbar(sm, ax=ax, location = 'left')
# cbar.set_label('\u03C3 [$N/mm^2$]', rotation=90)
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# plt.title("Axial Stress")
# ax.view_init(25 , -50)
# plt.xlabel("x")
# plt.ylabel("y")


# plt.savefig('Stress_Beam_Crane.pdf',format="pdf")
# plt.show()

# ###############################Deformation-x###################################
# fig5 = plt.figure()
# ax = fig5.add_subplot(111, projection='3d')

# Ux = U[:,0:1]
# Ux_mat = np.concatenate((Ux, np.zeros((n_nodes,5))),1)
# nodes_x_disp_coord_scaled = nodes_coord + scale*Ux_mat
# lin = 50 #linspace length

# i = 0
# for i in range(len(elements[0])):
#     if Ux[int(elements[1][i,0]-1)] > Ux[int(elements[1][i,1]-1)]:
#         Dmin = (Ux[int(elements[1][i,1]-1)]-min(Ux))/(max(Ux)-min(Ux)) #node with smallest displacement
#         Dmax = (Ux[int(elements[1][i,0]-1)]-min(Ux))/(max(Ux)-min(Ux)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes_x_disp_coord_scaled[int(elements[1][i,1]-1),0], nodes_x_disp_coord_scaled[int(elements[1][i,0]-1),0],lin)
#         ys = np.linspace(nodes_x_disp_coord_scaled[int(elements[1][i,1]-1),1], nodes_x_disp_coord_scaled[int(elements[1][i,0]-1),1],lin)
#         zs = np.linspace(nodes_x_disp_coord_scaled[int(elements[1][i,1]-1),2], nodes_x_disp_coord_scaled[int(elements[1][i,0]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

#     elif Ux[int(elements[1][i,0]-1)] <= Ux[int(elements[1][i,1]-1)]:
#         Dmin = (Ux[int(elements[1][i,0]-1)]-min(Ux))/(max(Ux)-min(Ux)) #node with smallest displacement
#         Dmax = (Ux[int(elements[1][i,1]-1)]-min(Ux))/(max(Ux)-min(Ux)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes_x_disp_coord_scaled[int(elements[1][i,0]-1),0], nodes_x_disp_coord_scaled[int(elements[1][i,1]-1),0],lin)
#         ys = np.linspace(nodes_x_disp_coord_scaled[int(elements[1][i,0]-1),1], nodes_x_disp_coord_scaled[int(elements[1][i,1]-1),1],lin)
#         zs = np.linspace(nodes_x_disp_coord_scaled[int(elements[1][i,0]-1),2], nodes_x_disp_coord_scaled[int(elements[1][i,1]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

# j = 0        
# for j in range(n_nodes-1):
#     ax.scatter(nodes_x_disp_coord_scaled[j,0],nodes_x_disp_coord_scaled[j,1],nodes_x_disp_coord_scaled[j,2],s=4, color = 'black')

# i = 0
# for i in range(n_nodes-1):
#     if (nodes[4][i]).any != 0:
#         coordF = nodes_x_disp_coord_scaled[i,:]
#         ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.2, color = 'black')

# norm = mpl.colors.Normalize(min(Ux), max(Ux))
# cmap = plt.get_cmap('jet', 100)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig5.colorbar(sm, ax=ax, location = 'left')
# cbar.set_label('$U_x$ [$mm$]', rotation=90)
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# plt.title("Displacement in x-direction")
# ax.view_init(25 , -50)
# plt.xlabel("x")
# plt.ylabel("y")


# plt.savefig('x_Displacement_Beam_Crane.pdf',format="pdf")
# plt.show()

# ###############################Deformation-y###################################
# fig6 = plt.figure()
# ax = fig6.add_subplot(111, projection='3d')

# Uy = U[:,1:2]
# Uy_mat = np.concatenate((np.zeros((len(nodes[0]),1)), Uy, np.zeros((n_nodes,4))),1)
# nodes_y_disp_coord_scaled = nodes_coord + scale*Uy_mat
# lin = 50 #linspace length

# i = 0
# for i in range(len(elements[0])):
#     if Uy[int(elements[1][i,0]-1)] > Uy[int(elements[1][i,1]-1)]:
#         Dmin = (Uy[int(elements[1][i,1]-1)]-min(Uy))/(max(Uy)-min(Uy)) #node with smallest displacement
#         Dmax = (Uy[int(elements[1][i,0]-1)]-min(Uy))/(max(Uy)-min(Uy)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes_y_disp_coord_scaled[int(elements[1][i,1]-1),0], nodes_y_disp_coord_scaled[int(elements[1][i,0]-1),0],lin)
#         ys = np.linspace(nodes_y_disp_coord_scaled[int(elements[1][i,1]-1),1], nodes_y_disp_coord_scaled[int(elements[1][i,0]-1),1],lin)
#         zs = np.linspace(nodes_y_disp_coord_scaled[int(elements[1][i,1]-1),2], nodes_y_disp_coord_scaled[int(elements[1][i,0]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

#     elif Uy[int(elements[1][i,0]-1)] <= Uy[int(elements[1][i,1]-1)]:
#         Dmin = (Uy[int(elements[1][i,0]-1)]-min(Uy))/(max(Uy)-min(Uy)) #node with smallest displacement
#         Dmax = (Uy[int(elements[1][i,1]-1)]-min(Uy))/(max(Uy)-min(Uy)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes_y_disp_coord_scaled[int(elements[1][i,0]-1),0], nodes_y_disp_coord_scaled[int(elements[1][i,1]-1),0],lin)
#         ys = np.linspace(nodes_y_disp_coord_scaled[int(elements[1][i,0]-1),1], nodes_y_disp_coord_scaled[int(elements[1][i,1]-1),1],lin)
#         zs = np.linspace(nodes_y_disp_coord_scaled[int(elements[1][i,0]-1),2], nodes_y_disp_coord_scaled[int(elements[1][i,1]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

# j = 0        
# for j in range(n_nodes-1):
#     ax.scatter(nodes_y_disp_coord_scaled[j,0],nodes_y_disp_coord_scaled[j,1],nodes_y_disp_coord_scaled[j,2],s=4, color = 'black')

# i = 0
# for i in range(n_nodes-1):
#     if (nodes[4][i]).any != 0:
#         coordF = nodes_y_disp_coord_scaled[i,:]
#         ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.2, color = 'black')

# norm = mpl.colors.Normalize(min(Uy), max(Uy))
# cmap = plt.get_cmap('jet', 100)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig6.colorbar(sm, ax=ax, location = 'left')
# cbar.set_label('$U_y$ [$mm$]', rotation=90)
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# plt.title("Displacement in y-direction")
# ax.view_init(25 , -50)
# plt.xlabel("x")
# plt.ylabel("y")


# plt.savefig('y_Displacement_Beam_Crane.pdf',format="pdf")
# plt.show()

# ###############################Deformation-z###################################
# fig7 = plt.figure()
# ax = fig7.add_subplot(111, projection='3d')

# Uz = -U[:,2:3]
# Uz_mat = np.concatenate((np.zeros((n_nodes,2)), Uz, np.zeros((n_nodes,3))),1)
# nodes_z_disp_coord_scaled = nodes_coord + scale*Uz_mat
# lin = 50 #linspace length

# i = 0
# for i in range(len(elements[0])):
#     if Uz[int(elements[1][i,0]-1)] > Uz[int(elements[1][i,1]-1)]:
#         Dmin = (Uz[int(elements[1][i,1]-1)]-min(Uz))/(max(Uz)-min(Uz)) #node with smallest displacement
#         Dmax = (Uz[int(elements[1][i,0]-1)]-min(Uz))/(max(Uz)-min(Uz)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes_z_disp_coord_scaled[int(elements[1][i,1]-1),0], nodes_z_disp_coord_scaled[int(elements[1][i,0]-1),0],lin)
#         ys = np.linspace(nodes_z_disp_coord_scaled[int(elements[1][i,1]-1),1], nodes_z_disp_coord_scaled[int(elements[1][i,0]-1),1],lin)
#         zs = np.linspace(nodes_z_disp_coord_scaled[int(elements[1][i,1]-1),2], nodes_z_disp_coord_scaled[int(elements[1][i,0]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

#     elif Uz[int(elements[1][i,0]-1)] <= Uz[int(elements[1][i,1]-1)]:
#         Dmin = (Uz[int(elements[1][i,0]-1)]-min(Uz))/(max(Uz)-min(Uz)) #node with smallest displacement
#         Dmax = (Uz[int(elements[1][i,1]-1)]-min(Uz))/(max(Uz)-min(Uz)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes_z_disp_coord_scaled[int(elements[1][i,0]-1),0], nodes_z_disp_coord_scaled[int(elements[1][i,1]-1),0],lin)
#         ys = np.linspace(nodes_z_disp_coord_scaled[int(elements[1][i,0]-1),1], nodes_z_disp_coord_scaled[int(elements[1][i,1]-1),1],lin)
#         zs = np.linspace(nodes_z_disp_coord_scaled[int(elements[1][i,0]-1),2], nodes_z_disp_coord_scaled[int(elements[1][i,1]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

# j = 0        
# for j in range(n_nodes-1):
#     ax.scatter(nodes_z_disp_coord_scaled[j,0],nodes_z_disp_coord_scaled[j,1],nodes_z_disp_coord_scaled[j,2],s=4, color = 'black')

# i = 0
# for i in range(n_nodes-1):
#     if (nodes[4][i]).any != 0:
#         coordF = nodes_z_disp_coord_scaled[i,:]
#         ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.2, color = 'black')

# norm = mpl.colors.Normalize(min(Uz),max(Uz))
# cmap = plt.get_cmap('jet', 100)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig7.colorbar(sm, ax=ax, location = 'left')
# cbar.set_label('$U_z$ [$mm$]', rotation=90)
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# plt.title("Absolute Displacement in z-direction")
# ax.view_init(25 , -50)
# plt.xlabel("x")
# plt.ylabel("y")


# plt.savefig('z_Displacement_Beam_Crane.pdf',format="pdf")
# plt.show()

# #############################Deformation-total#################################
# fig8 = plt.figure()
# ax = fig8.add_subplot(111, projection='3d')

# U_tot = np.zeros(len(U))
# i = 0
# for i in range(len(U)):
#     U_tot[i] = np.sqrt(U[i,0]**2+U[i,1]**2+U[i,2]**2)
# lin = 50 #linspace length

# i = 0
# for i in range(len(elements[0])):
#     if U_tot[int(elements[1][i,0]-1)] > U_tot[int(elements[1][i,1]-1)]:
#         Dmin = (U_tot[int(elements[1][i,1]-1)]-min(U_tot))/(max(U_tot)-min(U_tot)) #node with smallest displacement
#         Dmax = (U_tot[int(elements[1][i,0]-1)]-min(U_tot))/(max(U_tot)-min(U_tot)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes[2][int(elements[1][i,1]-1),0], nodes[2][int(elements[1][i,0]-1),0],lin)
#         ys = np.linspace(nodes[2][int(elements[1][i,1]-1),1], nodes[2][int(elements[1][i,0]-1),1],lin)
#         zs = np.linspace(nodes[2][int(elements[1][i,1]-1),2], nodes[2][int(elements[1][i,0]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

#     elif U_tot[int(elements[1][i,0]-1)] <= U_tot[int(elements[1][i,1]-1)]:
#         Dmin = (U_tot[int(elements[1][i,0]-1)]-min(U_tot))/(max(U_tot)-min(U_tot)) #node with smallest displacement
#         Dmax = (U_tot[int(elements[1][i,1]-1)]-min(U_tot))/(max(U_tot)-min(U_tot)) #node with largest displacement
#         grad = np.linspace(Dmin,Dmax,lin)
#         xs = np.linspace(nodes[2][int(elements[1][i,0]-1),0], nodes[2][int(elements[1][i,1]-1),0],lin)
#         ys = np.linspace(nodes[2][int(elements[1][i,0]-1),1], nodes[2][int(elements[1][i,1]-1),1],lin)
#         zs = np.linspace(nodes[2][int(elements[1][i,0]-1),2], nodes[2][int(elements[1][i,1]-1),2],lin)
        
#         j = 0
#         for j in range(lin-1):
#             ax.plot(xs[j:j+2], ys[j:j+2], zs[j:j+2], color=plt.cm.jet(grad[j]), linewidth = 1.0)

# j = 0        
# for j in range(n_nodes-1):
#     ax.scatter(nodes[2][j,0],nodes[2][j,1],nodes[2][j,2],s=4, color = 'black')

# i = 0
# for i in range(n_nodes-1):
#     if (nodes[4][i]).any != 0:
#         coordF = nodes[2][i,:]
#         ax.quiver(coordF[0],coordF[1],coordF[2],nodes[4][i,0],nodes[4][i,1],nodes[4][i,2], length=0.2, color = 'black')

# norm = mpl.colors.Normalize(min(U_tot), max(U_tot))
# cmap = plt.get_cmap('jet', 100)
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig8.colorbar(sm, ax=ax, location = 'left')
# cbar.set_label('$U_{total}$ [$mm$]', rotation=90)
# plt.xlim(-2000, 18000)
# plt.ylim(-2000, 18000)
# plt.title("Total Displacement")
# ax.view_init(25 , -50)
# plt.xlabel("x")
# plt.ylabel("y")


# plt.savefig('Total_Displacement_Beam_Crane.pdf',format="pdf")
# plt.show()



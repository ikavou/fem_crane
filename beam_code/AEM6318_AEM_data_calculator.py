#Ioannis Kavoukis, 6318
#AEM - data for crane problem

import numpy as np

def AEM_dataCalc(AEM):
    
    AEM = str(AEM)
    a = int(AEM[0])
    b = int(AEM[1])
    c = int(AEM[2])
    d = int(AEM[3])
    
    F = 20000 #N
    E = 210000 #N/mm2
    
    A = (1.2*(1+0.1*a+0.01*b))*1000 #mm
    B = (1+0.1*b+0.01*c)*1000 #mm
    phi = np.deg2rad(60+(b+0.1*d)) #rad
    theta = np.deg2rad(45+(a+0.1*c)) #rad
    A0=(6*(0.5+(10*d+c)/100))*100 #mm2
    A_diag = 0.5*A0
    A_bar = 1.5*A0
    L=(1.5*(1+0.1*d+0.01*c))*1000 #mm
    v = 0.3 #poisson
    if d == 0:
        n = 5
    else:
        n = 5 + np.mod(a,d)
    
    
    return [F,E,A,B,phi,theta,A_diag,A_bar,L,n,v]
    
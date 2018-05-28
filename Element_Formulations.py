import numpy as np
from FEAFunctions.FEAFunctions import FEA as fea
fea=fea()

'''
This script gives reaction force in one single hexahedral 3D 8 node element using three different element formulations:

- 8 integration points (full integration) using the principle of stationary potential energy. 
- 1 integration point (reduced integration) using the principle of stationary potential energy. 
- 8 integration points (full integration) using B-Bar method of the Hu Washizu three field formulation 
  with the modified B matrix formulated using the volume averaged volumetric strain to ease volumetric locking issues for nearly incompressible materials.

nodedisplacement specifies the displacement at the nodes. 

nodelocations specify the nodelocation of the element's nodes in numerical order according 
to a right handed x,y,z cartesian coordinate system
as demonstrated under. The nodedisplacement are specified in the same manner.

node:         |     1      |     2      |    3      |     4       |     5     |     6     |     7    |     8    |
coordinate:   | x   y   z  | x   y   z  | x  y   z  |  x  y   z   | x   y  z  | x   y  z  | x  y  z  |  x  y  z |
nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]

The integration points for the fully integrated formulations are numbered and located according to the closest node through the Gausspoints array.

Results are given in order of the nodes/integration points.

Strain and stress is defined in following order: 11, 22, 33, 12, 23, 31

'''

'''
Loadcases, comment in/out what to run
'''

# --------------- No volume change ---------------

# v=0.5 #poisson ratio
# dx=0.1
# dy=dz=-dx*v*0.5
# nodedisplacement=[0.0,-dy,-dz,  dx,-dy,-dz,  dx,dy,-dz,  0.0,dy,-dz,  0.0,-dy,dz, dx,-dy,dz, dx,dy,dz, 0.0,dy,dz] 
# nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]

# --------------- Hourglass ---------------

# nodedisplacement=[-dx,0.0,0.0,  -dx,0.0,0.0,  dx,0.0,0.0,  dx,0.0,0.0,  dx,0.0,0.0, dx,0.0,0.0, -dx,0.0,0.0, -dx,0.0,0.0]

# -------------- Pressurized pipe to benchmark volumetric locking ----------------

nodelocations=[3.25,0.,0.,  4.,0.,0.,  3.923141,0.780361,0.,  3.187552,0.634044,0.,  3.25,0.,1.,  4.,0.,1.,  3.923141,0.780361,1.,  3.187552,0.634044,1.]

defnodelocations=[ 3.26582, 0, 0,    4.01289, 0, 0,    3.93579, 0.782877, 0,     3.20306, 0.637129, 0,      3.26582, 0, 1,     4.01289, 0, 1,     3.93579, 0.782877, 1,     3.20306, 0.637129, 1 ] #C3D8 no

nodedisp=[]
for i in range(len(nodelocations)):
    n=-nodelocations[i]+defnodelocations[i]
    nodedisp.append(float(n))
    
nodedisplacement=nodedisp   

'''
Material Parameters
'''

v=0.45
E=100
K=E/(3*(1-(2*v)))
G=(3.*K*E)/((9.*K)-E)
C=fea.stiffnessmatrix(v=v,E=E)
Volume=fea.Volume(nodelocations=nodelocations)

'''
Element formulations, print what you want by modifying the formulations, comment in/out what to run.
'''

def main():
    OneFieldFullInt()
    OneFieldRedInt()
    ThreeFieldBBarLocking()

#----------- End input -------------------

def OneFieldFullInt():
    '''
    This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method.
    '''
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)] , [1./(3**0.5),-1./(3**0.5),-1./(3**0.5)]  ,   [1./(3**0.5),1./(3**0.5),-1./(3**0.5)],  [-1./(3**0.5),1./(3**0.5),-1./(3**0.5)], [-1./(3**0.5),-1./(3**0.5),1./(3**0.5)], [1./(3**0.5),-1./(3**0.5),1./(3**0.5)],   [1./(3**0.5),1./(3**0.5),1./(3**0.5)], [-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
    Ke=np.zeros((24,24))
    S=[]
    e=[]
    for Gauss in Gausspoints:
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        Bt=B.transpose()
        Ke=(w*Bt*C*B*detJ)+Ke
        eu=B*(np.matrix(nodedisplacement).transpose())
        e.append(eu)
        s=C*eu
        S.append(s)
    F=Ke*(np.matrix(nodedisplacement).transpose())
    print("--------------------------- Standard 8 node element formulation with eight integration points (fully integrated) ----------------------------------")
    print("Reaction forces:")
    print(F)
    print("Stress:")
    print(S)
    print("Strain:")
    print(e)

    # --------------- Calculate rank defiency ---------------
    
    orderKe=len(Ke)
    nRB=6
    rankKe=np.linalg.matrix_rank(Ke)
    properrank=orderKe-nRB
    rankdef=float(properrank)-float(rankKe)
    print("Rank defiancy: ", rankdef)

def OneFieldRedInt():

    '''
    This section calculates reaction force in a reduced integration 8 node hexahedron with 1 integration point using the standard Gaussian integration method
    '''
    
    Gausspoints=np.array([[0.,0.,0.]])
    w=8. #Sum of all weights must be 8. ref colorado net course.
    Ke=np.zeros((24,24))
    S=[]
    e=[]
    for Gauss in Gausspoints:
        
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        eu=B*(np.matrix(nodedisplacement).transpose())
        s=C*eu
        S.append(s)  
        e.append(eu)
        Bt=B.transpose()
        Ke=(w*Bt*C*B*detJ)+Ke
    F=Ke*(np.matrix(nodedisplacement).transpose())
    print("--------------------------- Standard eight node element formulation with one integration point (reduced integration) ----------------------------------")
    print("Reaction forces:")
    print(F)
    print("Stress:")
    print(S)
    print("Strain:")
    print(e)   
    
    # --------------- Calculate rank defiency ---------------
    
    orderKe=len(Ke)
    nRB=6
    rankKe=np.linalg.matrix_rank(Ke)
    properrank=orderKe-nRB
    rankdef=float(properrank)-float(rankKe)
    print("Rank defiancy: ", rankdef)

def ThreeFieldBBarLocking():

    '''
    This section calculates reaction force, stress and strain in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the B-Bar method of the Hu Washizu three field formulation with diff between volumetric strain in point and avg vol strain plus the actual strain in point 
    as added degrees of freedom to illustrate the locking effect and how the B-Bar can handle this.
    '''
    Bv=np.zeros((6,24))
    Be=np.zeros((6,24))
    Ke=np.zeros((24,24))      
    
    e=[]
    s=[]
    m=(np.matrix([1.,1.,1.,0.,0.,0.])).transpose()
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],
                          [-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],
                          [1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
    for Gauss in Gausspoints:

        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])        
        Bv=((m*m.transpose()*B*w*detJ))+Bv              #B-matrix giving volume averaged (therefore the integration scheme) volumetric strain.             
                                                          
    Bv=Bv*(1./Volume)                                   #Divide by total volume to get average.              
    ev=Bv*(np.matrix(nodedisplacement).transpose())     #Volumetric strain
    
    for Gauss in Gausspoints:
        
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])        
        Be=(B+((1./3.)*(Bv-(m*m.transpose()*B))))                     #Modified B matrix                    

        eu=Be*(np.matrix(nodedisplacement).transpose())               #Assumed strain in integration point              
        e.append(eu)                                                  #Append assumed strains in integration points             
        
        S=C*eu                                                        #Assumed stress in integration point  
        s.append(S)                                                   #Append assumed stress integration points       
        
        Ke=(Be.transpose()*C*Be*w*detJ)+Ke                            #Assemble element stiffness matrix, Ke.
    Fs=Ke*(np.matrix(nodedisplacement).transpose())                   #Get reaction forces  
    
    print("--------------------------- B-bar method on eight node element with eight integration points (fully integrated) ----------------------------------")
    print("Reaction forces:")
    print(Fs)
    print("B-Bar Stress field:")
    print(s)
    print("B-Bar strain field:")
    print(e)
    print("Volumetric strain:")
    print(ev[0,0])

main()
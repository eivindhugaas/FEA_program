import numpy as np

from FEAFunctions.FEAFunctions import FEA as fea
fea=fea()

dx=0.1 #mm
dy=-0.025#-0.015 #mm 
dz=-0.025#-0.015 #mm 
 
nodedisplacement=[0.0,-dy,-dz,  dx,-dy,-dz,  dx,dy,-dz,  0.0,dy,-dz,  0.0,-dy,dz, dx,-dy,dz, dx,dy,dz, 0.0,dy,dz] #displacement along dx with free planes on y and z dir. This is u with hat = u on Su, where Su is the surfaces whichj displace. This part is normally done by connectivity matrixes.
#nodedisplacement=[-dx,0.0,0.0,  -dx,0.0,0.0,  dx,0.0,0.0,  dx,0.0,0.0,  dx,0.0,0.0, dx,0.0,0.0, -dx,0.0,0.0, -dx,0.0,0.0] # Hourglass in x direction.
d=1-0.99665551
nodedisplacement=[d,d,d,  -d,d,d,  -d,-d,d, d,-d,d,  d,d,-d,  -d,d,-d,  -d,-d,-d, d,-d,-d] #1% reduction of volume
nodelocations=[-1.,-23.,-1.,  1.,-23.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-2.,1.,  1.,-2.,1.,  1.,1.,1.,  -1.,1.,1.]
#nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
v=0.3
E=210000. #MPa
K=160000. #Mpa A material with a bulk modulus of 35000 MPa loses one percent of its volume when subjected to an external pressure of 0.35 GPa = 350MPa (~3500 bar).

C=fea.stiffnessmatrix(v=v,E=E)

#---------- End input ----------#

'''
This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method
'''
Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
w=1.

Ke=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
S=[]
for Gauss in Gausspoints:
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    eu=B*(np.matrix(nodedisplacement).transpose())
    s=C*eu
    S.append(s)
    Bt=B.transpose()
    Ke=(w*Bt*C*B*detJ)+Ke
F=Ke*(np.matrix(nodedisplacement).transpose())
Fr=F
print("--------------------------- eight integration points ----------------------------------")
print("Reaction forces fully integrated:")
print(F[0,0])
print("Stress:")
print(S[0][0])

#calculate rank defiency

orderKe=len(Ke)
nRB=6
rankKe=np.linalg.matrix_rank(Ke)
properrank=orderKe-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)

'''
This section calculates reaction force in a reduced integration 8 node hexahedron with 1 integration point using the standard Gaussian integration method
'''

Gausspoints=np.array([[0.,0.,0.]])
w=8. #Sum of all weights must be 8. ref colorado net course.

Ke=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Stress=[]
for Gauss in Gausspoints:
    
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    eu=B*(np.matrix(nodedisplacement).transpose())
    s=C*eu
    Stress.append(s)    
    Bt=B.transpose()
    Ke=(w*Bt*C*B*detJ)+Ke
F=Ke*(np.matrix(nodedisplacement).transpose())
print("--------------------------- one integration point ----------------------------------")
print("Reaction forces fully integrated:")
print(F[0,0])
print("Stress:")
print(Stress[0][0])
#calculate rank defiency

orderKe=len(Ke)
nRB=6
rankKe=np.linalg.matrix_rank(Ke)
properrank=orderKe-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)

'''
This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
and the Hellinger Reissner two field formulation with hydrostatic stress as added degree of freedom. Takes into account that the material is compressible with a K value.
'''
Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
w=1.

Ke=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Kus=[[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]]
Kss=0.
Stress=[]
Stressdev=[]
m=(np.matrix([1.,1.,1.,0.,0.,0.])).transpose()
I0=np.matrix(([1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0.5,0,0],[0,0,0,0,0.5,0],[0,0,0,0,0,0.5]))
#Y=np.dot(m,(m.transpose))
Cd=(2*K*(I0-((1/3)*m*(m.transpose()))))
for Gauss in Gausspoints:
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Bt=B.transpose()
    Ke=(w*Bt*Cd*B*detJ)+Ke

    Ns=fea.HydrostatStresshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Kus=(Bt*m*Ns*w*detJ)+Kus
    Kss=(Ns*Ns*(1/K)*w*detJ)+Kss
    
K=Ke+(Kus*(1./Kss)*(Kus.transpose()))
#K=Ke+(Kus*(Kss)*(Kus).transpose())
Sh=float((Kus).transpose()*(1/Kss)*(np.matrix(nodedisplacement).transpose()))

Fs=K*(np.matrix(nodedisplacement).transpose())
    
print("--------------------------- two field, hydrostatic stress as second degree of freedom ----------------------------------")
print("Reaction forces fully integrated:")
print(Fs[0,0])
print("Stress:")
print(Sh)

#calculate rank defiency

orderKe=len(K)
nRB=6
rankKe=np.linalg.matrix_rank(K)
properrank=orderKe-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)


'''
This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
and the Hellinger Reissner two field formulation with average stress as added degree of freedom. Due to average stress being the degree of freedom, there are less stress DOF than U-DOF.
This should not be allowed, however, if viewed as each node has the average stress as degree of freedom with 1 as interpolation function, there are in theory twice as many stress DOF as U DOF.
'''
Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
w=1.

Ksu=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Kss=[[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.]]
Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
s1=Ns.shape
B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
s2=B.shape
size=(s1[1],s2[1])
Ksu=np.zeros(size)

size=(48,48)
Kes=np.zeros(size)

size=(48,48)
Kee=np.zeros(size)

size=(s1[1],s1[1])
Kss=np.zeros(size)

K=np.zeros((24,24))

for Gauss in Gausspoints:
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Bt=B.transpose()

    Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Ksu=(Ns.transpose()*B*w*detJ)+Ksu
    Kss=-(Ns.transpose()*np.linalg.inv(C)*Ns*w*detJ)+Kss
    
K=((Ksu.transpose())*(np.linalg.inv(Kss))*Ksu)  #Normal Ke
Fs=-K*(np.matrix(nodedisplacement).transpose()) #Normal Reaction force
S=-np.linalg.inv((Kss.transpose()))*Ksu*(np.matrix(nodedisplacement).transpose())
print("--------------------------- two field, average stress as second degree of freedom ----------------------------------")
print("Reaction forces fully integrated:")
print(Fs[0,0])
print("Stress:")
print(S[0,0])

#calculate rank defiency

orderKe=len(K)
nRB=6
rankKe=np.linalg.matrix_rank(K)
properrank=orderKe-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)


'''
This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
and the Hu Washizu three field formulation with Average stress and average strain as added degree of freedom.
'''

Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
w=1.

Ke=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Ksu=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Kee=[[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.]]
Kes=[[0.],[0.],[0.],[0.],[0.],[0.]]
Kss=0.



Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
s1=Ns.shape
B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
s2=B.shape
size=(s1[0],s2[1])
Ksu=np.zeros(size)


size=s1
Kes=np.zeros(size)

size=s1
Kee=np.zeros(size)

Stress=[]
Stressdev=[]
m=(np.matrix([1,1,1,0,0,0])).transpose()

for Gauss in Gausspoints:
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Bt=B.transpose()
    Ke=(w*Bt*C*B*detJ)+Ke
    

    Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Ne=fea.AvgStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    
    Ksu=((Ns.transpose())*B*w*detJ)+Ksu   # gives strain from disp DOF
    Kes=-((Ne.transpose())*Ns*w*detJ)+Kes   # gives stress in strain shape from stress DOF
    Kee=(Ne.transpose()*C*Ne*w*detJ)+Kee  # gives stress in strain shape from strain DOF
    
K=(Ksu.transpose())*(np.linalg.inv(Kes))*Kee*((np.linalg.inv(Kes).transpose()))*Ksu
e=-((np.linalg.inv(Kes).transpose()))*Ksu*(np.matrix(nodedisplacement).transpose())
S=-((np.linalg.inv(Kes).transpose()))*Kee*e


Fs=K*(np.matrix(nodedisplacement).transpose())
    
print("--------------------------- three field, average stress as second degree of freedom and average strain as third ----------------------------------")
print("Reaction forces fully integrated:")
print(Fs[0,0])
print("Stress:")
print(S[0,0])
print("Strain:")
print(e[0,0])

#calculate rank defiency

orderKe=len(K)
nRB=6
rankKe=np.linalg.matrix_rank(K)
properrank=orderKe-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)


'''
This section calculates reaction force, stress and strain in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
and the Hu Washizu three field formulation with conventionally interpolated stress and strain as added degrees of freedom.
'''

Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
w=1.

Ke=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Ksu=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Kee=[[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.]]
Kes=[[0.],[0.],[0.],[0.],[0.],[0.]]

#Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
Ns=Ne=fea.norderStrainshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
Ns=Ne

s1=Ns.shape
B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
s2=B.shape
size=(s1[1],s2[1])

Ksu=np.zeros(size)
Ne.shape
s3=Ns.shape

size=(s1[1],s1[1])
Kes=np.zeros(size)

size=(s1[1],s1[1])
Kee=np.zeros(size)


for Gauss in Gausspoints:
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Bt=B.transpose()
    Neexp=fea.twoorderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Nexp=fea.norderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Ne=Neexp
    Ns=Nexp  
    Ksu=((Ns.transpose())*B*w*detJ)+Ksu   # gives strain from disp DOF
    Kes=-((Ne.transpose())*Ns*w*detJ)+Kes   # gives stress in strain shape from stress DOF
    Kee=(Ne.transpose()*C*Ne*w*detJ)+Kee  # gives stress in strain shape from strain DOF


K=(Ksu.transpose())*(np.linalg.inv(Kes))*Kee*((np.linalg.inv(Kes).transpose()))*Ksu
e=-((np.linalg.inv(Kes).transpose()))*Ksu*(np.matrix(nodedisplacement).transpose())
S=-((np.linalg.inv(Kes).transpose()))*Kee*e

Fs=K*(np.matrix(nodedisplacement).transpose())
 
    
print("--------------------------- three field, conventinonally interpolated stress as second degree of freedom and 2nd order equation interpolated strain as third ----------------------------------")
print("Reaction forces fully integrated:")
print(Fs[0,0])
print("Stress:")
print(S[0,0])
print("Strain:")
print(e[0,0])
#calculate rank defiency

orderK=len(K)
nRB=6
rankKe=np.linalg.matrix_rank(K)
properrank=orderK-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)

'''
This section calculates reaction force, stress and strain in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
and the B-Bar method of the Hu Washizu three field formulation with conventionally interpolated stress and strain as added degrees of freedom.
'''

Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
w=1.

Ke=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Ksu=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
Kee=[[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.]]
Kes=[[0.],[0.],[0.],[0.],[0.],[0.]]

Ns=Ne=fea.AvgStrainshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
Ns=Ne
#Ns=fea.norderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
s1=Ns.shape
B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
s2=B.shape
size=(s1[1],s2[1])
s4=Ne.shape
#size=(48,24)
Ksu=np.zeros(size)

Be=np.zeros((s2))
Ne.shape

s3=Ns.shape

size=(s1[0],s1[1])
Kes=np.zeros(size)

size=(s4[1],s4[1])
Kee=np.zeros(size)


Kss=0.
Stress=[]
Stressdev=[]
m=(np.matrix([1,1,1,0,0,0])).transpose()
A=fea.UStraintoAvgStrainshapefunc()
eavg=np.matrix([0.,0.,0.,0.,0.,0.]).transpose()
for Gauss in Gausspoints:
    
    #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Be=(A*B)+Be#flips the B matrix into giving an addition as average to the avg. strain.
    
    Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Neexp=fea.norderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    Ne=fea.AvgStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    #Ns=Neexp
   
    Ksu=((Ns.transpose())*B*w*detJ)+Ksu   # gives strain from disp DOF
    Kes=-((Ne.transpose())*Ns*w*detJ)+Kes   # gives stress in strain shape from stress DOF
    Kee=(Ne.transpose()*C*Ne*w*detJ)+Kee  # gives stress in strain shape from strain DOF
    


eavg=Be*(np.matrix(nodedisplacement).transpose())

K=-(Ksu.transpose())*(np.linalg.inv(Kes))*Kee*Be*detJ
e=Be*(np.matrix(nodedisplacement).transpose())
S=-(np.linalg.inv(Kes))*Kee*Be*(np.matrix(nodedisplacement).transpose())

Fs=K*(np.matrix(nodedisplacement).transpose())
    
print("--------------------------- three field, conventinonally interpolated stress as second degree of freedom and conventionally interpolated strain as third ----------------------------------")
print("Reaction forces fully integrated:")
print(Fs[0,0])
print("Stress:")
print(S[0,0])
print("Strain:")
print(e[0,0])
#calculate rank defiency

orderK=len(K)
nRB=6
rankKe=np.linalg.matrix_rank(K)
properrank=orderK-nRB
rankdef=float(properrank)-float(rankKe)
print("Rank defiancy: ", rankdef)


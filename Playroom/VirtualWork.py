import numpy as np

from FEAFunctions.FEAFunctions import FEA as fea
fea=fea()
v=0.5
dx=0.1 #mm
dy=-0.025
dz=-0.025
dy=dz=-dx*v*0.5
nodedisplacement=[0.0,-dy,-dz,  dx,-dy,-dz,  dx,dy,-dz,  0.0,dy,-dz,  0.0,-dy,dz, dx,-dy,dz, dx,dy,dz, 0.0,dy,dz] #displacement along dx with free planes on y and z dir. This is u with hat = u on Su, where Su is the surfaces whichj displace. This part is normally done by connectivity matrixes.
v=0.499999999

nodedisplacement=[-dx,0.0,0.0,  -dx,0.0,0.0,  dx,0.0,0.0,  dx,0.0,0.0,  dx,0.0,0.0, dx,0.0,0.0, -dx,0.0,0.0, -dx,0.0,0.0] # Hourglass in x direction.

nodelocations=[3.25,0.,0.,  4.,0.,0.,  3.923141,0.780361,0.,  3.187552,0.634044,0.,  3.25,0.,1.,  4.,0.,1.,  3.923141,0.780361,1.,  3.187552,0.634044,1.]

#defnodelocations=[3.41301, 0, 0,   4.13293, 0, 0,    4.05352, 0.806296, 0,   3.34743, 0.665846, 0,     3.41301, 0, 1,          4.13293, 0, 1,                4.05352, 0.806296, 1,       3.34743, 0.665846, 1] #C3D20R
defnodelocations= [3.40816, 0, 0,   4.12894, 0, 0,    4.0496, 0.805516, 0,    3.34267, 0.664899, 0,     3.40816, 0, 1,          4.12894, 0, 1,                4.0496, 0.805516,1,         3.34267, 0.664899, 1] #C3D8
#defnodelocations= [3.40816+10, 0, 0,   4.12894+10, 0, 0,    4.0496+10, 0.805516, 0,    3.34267+10, 0.664899, 0,     3.40816+10, 0, 1,          4.12894+10, 0, 1,                4.0496+10, 0.805516,1,         3.34267+10, 0.664899, 1] 
#nodelocations=[-0.3401725, -0.35360125, -0.5, 0.4098275, -0.35360125, -0.5, 0.3329675000000001, 0.42675975, -0.5, -0.4026225000000001, 0.2804427500000001, -0.5, -0.3401725, -0.35360125, 0.5, 0.4098275, -0.35360125, 0.5, 0.3329675000000001, 0.42675975, 0.5, -0.4026225000000001, 0.2804427500000001, 0.5]
#defnodelocations=[-0.3047925, -0.363746875, -0.5, 0.41598749999999995, -0.363746875, -0.5, 0.33664749999999977, 0.441769125, -0.5, -0.37028250000000007, 0.301152125, -0.5, -0.3047925, -0.363746875, 0.5, 0.41598749999999995, -0.363746875, 0.5, 0.33664749999999977, 0.441769125, 0.5, -0.5254025000000002, 0.270297125, 0.5]
defnodelocations=[ 3.26582, 0, 0,    4.01289, 0, 0,    3.93579, 0.782877, 0,     3.20306, 0.637129, 0,      3.26582, 0, 1,     4.01289, 0, 1,     3.93579, 0.782877, 1,     3.20306, 0.637129, 1 ] #C3D8 no

nodedisp=[]
for i in range(len(nodelocations)):
    n=-nodelocations[i]+defnodelocations[i]
    nodedisp.append(float(n))
    
nodedisplacement=nodedisp   

print(nodedisp)

v=0.495
E=100 #MPa
K=E/(3*(1-(2*v)))
G=(3.*K*E)/((9.*K)-E)
C=fea.stiffnessmatrix(v=v,E=E)
Volume=fea.Volume(nodelocations=nodelocations)


#---------- End input ----------#

def main():
    OneFieldFullInt()
    #OneFieldRedInt()
    TwoFieldAvgStress()
    #TwoFieldLocking()
    #ThreeFieldAvgStressAvgStrain()
    #ThreeFieldBBarAvgStressAvgStrain()
    ThreeFieldBBarLocking()
    #ThreeFieldNormalStressThirdStrain()


def OneFieldFullInt():
    '''
    This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method
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
    Fr=F
    print("--------------------------- Standard element formulation with eight integration points ----------------------------------")
    print("Reaction forces:")
    print(F)

    print("Stress:")
    print(S)
    print("Strain:")
    print(e)
    #print(float(e[0][0]))   

    #calculate rank defiency
    
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
        
        #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        eu=B*(np.matrix(nodedisplacement).transpose())
        s=C*eu
        S.append(s)  
        e.append(eu)
        Bt=B.transpose()
        Ke=(w*Bt*C*B*detJ)+Ke
    F=Ke*(np.matrix(nodedisplacement).transpose())
    print("--------------------------- Standard element formulation with one integration point ----------------------------------")
    print("Reaction forces:")
    print(F)
    print("Stress:")
    print(S)
    print("Strain:")
    print(e)   
    #calculate rank defiency
    
    orderKe=len(Ke)
    nRB=6
    rankKe=np.linalg.matrix_rank(Ke)
    properrank=orderKe-nRB
    rankdef=float(properrank)-float(rankKe)
    print("Rank defiancy: ", rankdef)

def TwoFieldLocking():
    '''
    This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the Hellinger Reissner two field formulation with hydrostatic stress as added degree of freedom. Takes into account that the material is compressible with a K value.
    '''
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
    
    Kuu=np.zeros((24,24))
    Kus=np.zeros((24,48))
    Kss=0.
    m=(np.matrix([1.,1.,1.,0.,0.,0.])).transpose()
    I0=np.matrix(([1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0.5,0,0],[0,0,0,0,0.5,0],[0,0,0,0,0,0.5]))
    Cd=(2*G*(I0-((1/3)*m*(m.transpose()))))
    for Gauss in Gausspoints:
        #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        Bt=B.transpose()
        
        Kuu=(w*Bt*Cd*B*detJ)+Kuu
        
        Ns=fea.HydrostatStresshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        Kus=(Bt*m*Ns*w*detJ)+Kus
        
  
        
        Kss=(Ns.transpose()*Ns*(1/K)*w*detJ)+Kss
    
    print(Kss)
        
    Ke=Kuu+(Kus*(np.linalg.inv(Kss))*(Kus.transpose()))
    
    print(Ke)
    
    Fs=Ke*(np.matrix(nodedisplacement).transpose())
    
    print("--------------------------- Two field with hydrostatic stress as second degree of freedom ----------------------------------")
    print("Reaction forces:")
    print(Fs)
    print("Hydrostatic stress:")
   # print(float(Sh[0][0]))
   # print("Volumetric strain:")
   # print(float(ev[0][0]))        
        

def TwoFieldAvgStress():
    '''
    This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the Hellinger Reissner two field formulation with average stress as added degree of freedom. Due to average stress being the degree of freedom, there are less stress DOF than U-DOF.
    This should not be allowed, however, if viewed as each node has the average stress as degree of freedom with 1 as interpolation function, there are in theory twice as many stress DOF as U DOF.
    '''
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
    
    
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
    
    Ke=np.zeros((24,24))
    
    for Gauss in Gausspoints:
        #nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        Bt=B.transpose()
    
        Ns=fea.AvgStresshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        Ksu=(Ns.transpose()*B*w*detJ)+Ksu
        Kss=-(Ns.transpose()*np.linalg.inv(C)*Ns*w*detJ)+Kss

    Ke=((Ksu.transpose())*(np.linalg.inv(Kss))*Ksu)  #Normal Ke
    Fs=-Ke*(np.matrix(nodedisplacement).transpose()) #Normal Reaction force
    S=-np.linalg.inv((Kss.transpose()))*Ksu*(np.matrix(nodedisplacement).transpose())
    print("--------------------------- Two field with average stress as second degree of freedom ----------------------------------")
    print("Reaction forces fully integrated:")
    print(Fs)
    print("Average stress:")
    print(S[0,0])
    
def ThreeFieldAvgStressAvgStrain():
    '''
    This section calculates reaction force in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the Hu Washizu three field formulation with Average stress and average strain as added degree of freedom.
    '''
    
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.

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
    
    Ke=np.zeros((24,24))
    
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
        
    Ke=(Ksu.transpose())*(np.linalg.inv(Kes))*Kee*((np.linalg.inv(Kes).transpose()))*Ksu
    e=-((np.linalg.inv(Kes).transpose()))*Ksu*(np.matrix(nodedisplacement).transpose())
    S=-((np.linalg.inv(Kes).transpose()))*Kee*e
    
    
    Fs=Ke*(np.matrix(nodedisplacement).transpose())
        
    print("--------------------------- three field, average stress as second degree of freedom and average strain as third ----------------------------------")
    print("Reaction forces fully integrated:")
    print(Fs[0,0])
    print("Average stress:")
    print(S[0,0])
    print("Average strain:")
    print(e)


def ThreeFieldNormalStressThirdStrain():

    '''
    This section calculates reaction force, stress and strain in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the Hu Washizu three field formulation with conventionally interpolated stress and third order interpolated strain as added degrees of freedom.
    '''
    
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
       
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
        Neexp=fea.twoorderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2]) #threeorder
        Nexp=fea.norderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])    #normalorder
        Ne=Nexp #three order
        Ns=Nexp  #normal order
        Ksu=((Ns.transpose())*B*w*detJ)+Ksu   # gives strain from disp DOF
        Kes=-((Ne.transpose())*Ns*w*detJ)+Kes   # gives stress in strain shape from stress DOF
        Kee=(Ne.transpose()*C*Ne*w*detJ)+Kee  # gives stress in strain shape from strain DOF
    
    Ke=(Ksu.transpose())*(np.linalg.inv(Kes))*Kee*((np.linalg.inv(Kes).transpose()))*Ksu
    e=-((np.linalg.inv(Kes).transpose()))*Ksu*(np.matrix(nodedisplacement).transpose())
    S=-((np.linalg.inv(Kes).transpose()))*Kee*e
    
    Fs=Ke*(np.matrix(nodedisplacement).transpose())
     
    print("--------------------------- three field, conventinonally interpolated stress as second degree of freedom and 3d order equation interpolated strain as third ----------------------------------")
    print("Reaction forces fully integrated:")
    print(Fs[0,0])
    print("Stress:")
    print(S[0,0])
    print("Strain:")
    print(e)#[0,0])

def ThreeFieldBBarAvgStressAvgStrain():

    '''
    This section calculates reaction force, stress and strain in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the B-Bar method of the Hu Washizu three field formulation with avg interpolated stress and strain as added degrees of freedom.
    '''
    
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
    
    Ns=Ne=fea.AvgStrainshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
    Ns=Ne
    #Ns=fea.norderStrainshapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
    s1=Ns.shape
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
    
    size=(24,24)
    Ke=np.zeros(size)

    size=(6,24)
    Be=np.zeros(size)
    
    m=(np.matrix([1,1,1,0,0,0])).transpose()
    A=fea.UStraintoAvgStrainshapefunc()
    eavg=np.matrix([0.,0.,0.,0.,0.,0.]).transpose()
    
    for Gauss in Gausspoints:
        
        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])
        Be=(A*B)+Be#flips the B matrix into giving an addition as average to the avg. strain.
    
    for Gauss in Gausspoints:
    
        Ke=Be.transpose()*C*Be*w*detJ+Ke
        
    Fs=Ke*(np.matrix(nodedisplacement).transpose())
    
    e=Be*(np.matrix(nodedisplacement).transpose())

    S=C*e
    
    print("--------------------------- three field B-bar, average stress as second degree of freedom and average strain deducted using an enhanced B matrix from the displacement as third ----------------------------------")
    print("Reaction forces fully integrated:")
    print(Fs[0,0])
    print("Average stress:")
    print(float(S[0][0]))
    print("Average strain:")
    print(float(e[0][0]))


def ThreeFieldBBarLocking():

    '''
    This section calculates reaction force, stress and strain in a fully integrated 8 node hexahedron with 8 integration points using the standard Gaussian integration method 
    and the B-Bar method of the Hu Washizu three field formulation with diff between volumetric strain in point and avg vol strain plus the actual strain in point 
    as added degrees of freedom to illustrate the locking effect and how the B-Bar can handle this.
    '''
    
    Gausspoints=np.array([[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),-1./(3**0.5)],[1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),-1./(3**0.5)],[-1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),-1./(3**0.5),1./(3**0.5)],[1./(3**0.5),1./(3**0.5),1./(3**0.5)],[-1./(3**0.5),1./(3**0.5),1./(3**0.5)]])
    w=1.
    
    
    Ns=Ne=fea.AvgStrainshapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
    Ns=Ne
    s1=Ns.shape
    B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=0.,eta=0.,zeta=0.)
    s2=B.shape
    size=(s1[1],s2[1])
    s4=Ne.shape
    size=(6,24)
    
    Bv=np.zeros((s2))
    
    Be=np.zeros((s2))
    
    
    size=(24,24)
    Ke=np.zeros(size)      
    
    e=[]
    s=[]
    m=(np.matrix([1,1,1,0,0,0])).transpose()
    A=fea.UStraintoAvgStrainshapefunc()
    eavg=np.matrix([0.,0.,0.,0.,0.,0.]).transpose()
    m=(np.matrix([1.,1.,1.,0.,0.,0.])).transpose()
    I0=np.matrix(([1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0.5,0,0],[0,0,0,0,0.5,0],[0,0,0,0,0,0.5]))
    Cd=(2*G*(I0-((1/3)*m*(m.transpose()))))    
    
    for Gauss in Gausspoints:

        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])        
        
        Bv=((1./Volume)*(m*m.transpose()*B*w*detJ))+Bv  
    
    ev=Bv*(np.matrix(nodedisplacement).transpose())
    
    print("hello",ev)
    
    for Gauss in Gausspoints:

        B,detJ=fea.shapefunc(nodelocations=nodelocations,exi=Gauss[0],eta=Gauss[1],zeta=Gauss[2])        

        Be=(B+((1./3.)*(Bv-(m*m.transpose()*B))))                     #correct, gives diff between volumetric strain in point and avg vol strain plus the actual strain in point. So if ev=evavg Be gives e.
        
        eu=Be*(np.matrix(nodedisplacement).transpose())               #assumed strain field, damn smood.
        
        e.append(eu)                                                  #assumed strain field, damn smood and appended.
        
        S=(Cd*B*(np.matrix(nodedisplacement).transpose()))+((K/1)*ev) #Assumed stress field.
        
        s.append(S)                                                   #Assumed stress field appended.
        
        Ke=(Be.transpose()*C*Be*w*detJ)+Ke                        
    Fs=Ke*(np.matrix(nodedisplacement).transpose())
    
    print("--------------------------- three field B-bar, difference between strain and  strain deducted using an enhanced B matrix from the displacement as third ----------------------------------")
    print("Reaction forces fully integrated:")
    print(Fs)

    print("Stress, abstract:")
    print(float(s[0][0]))
    print("Strain, abstract:")
    print(float(e[0][0]))
    

main()
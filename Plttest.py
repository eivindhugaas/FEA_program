
import numpy as np
from sympy import *
import matplotlib as plt
'''
isoparametric, hexahedral, natural coordinates:
xi is an e with a krussedull on top
eta is an n with a long leg
zeta is a snake that is bowing, it is zzzzzzeeeetaaaaaa.
'''
Gauss=[-1./(3**0.5),-1./(3**0.5),-1./(3**0.5)]
Gauss=[1,1,-1]
nodelocations=[-1.,-1.,-1.,  1.,-1.,-1.,  1.,1.,-1.,  -1.,1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,  1.,1.,1.,  -1.,1.,1.]
exi=Gauss[0]
eta=Gauss[1]
zeta=Gauss[2]


t=0
nodal=[]
nodecoords=[]
for i in range(len(nodelocations)):
    n=i-t

    if n<2:

        nodal.append(nodelocations[i])
        
    else:
        nodal.append(nodelocations[i])
        nodecoords.append(nodal)
        nodal=[]
        t=i+1
    
#print(nodecoords)

#nodecoords=[[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.],[1.,-1.,1.],[1.,1.,1.],[-1.,1.,1.]]
ShapeFuncVectorx=[]
ShapeFuncVectory=[]
ShapeFuncVectorz=[]
ShapeFuncVectorxx=[]
ShapeFuncVectoryy=[]
ShapeFuncVectorzz=[]        
ShapeFuncVector=[]
for i in range(0,6,1):
    for coord in nodecoords:
        x=coord[0]
        y=coord[1]
        z=coord[2]
        #N=(1/8.)*(1.+(exi*x*x))*(1.+(eta*y*y))*(1.+(zeta*z*z))*((exi*x*x)+(eta*y*y)+(zeta*z*z)-2.)
        N=(1/8.)*(1.+(exi*x))*(1.+(eta*y))*(1.+(zeta*z))*((exi*x)+(eta*y)+(zeta*z)-2.)
        if i==0:
            ShapeFuncVectorx.append(N)
            ShapeFuncVectorx.append(0.)
            ShapeFuncVectorx.append(0.)
            ShapeFuncVectorx.append(0.)
            ShapeFuncVectorx.append(0.)
            ShapeFuncVectorx.append(0.)
            
        if i==1:
            ShapeFuncVectory.append(0.)                    
            ShapeFuncVectory.append(N)
            ShapeFuncVectory.append(0.)
            ShapeFuncVectory.append(0.)
            ShapeFuncVectory.append(0.)
            ShapeFuncVectory.append(0.)
        if i==2:
            ShapeFuncVectorz.append(0.)                    
            ShapeFuncVectorz.append(0.)
            ShapeFuncVectorz.append(N)
            ShapeFuncVectorz.append(0.)
            ShapeFuncVectorz.append(0.)
            ShapeFuncVectorz.append(0.)
        if i==3:
            ShapeFuncVectorxx.append(0.)
            ShapeFuncVectorxx.append(0.)
            ShapeFuncVectorxx.append(0.)
            ShapeFuncVectorxx.append(N)
            ShapeFuncVectorxx.append(0.)
            ShapeFuncVectorxx.append(0.)
        if i==4:
            ShapeFuncVectoryy.append(0.)
            ShapeFuncVectoryy.append(0.)
            ShapeFuncVectoryy.append(0.)
            ShapeFuncVectoryy.append(0.)                    
            ShapeFuncVectoryy.append(N)
            ShapeFuncVectoryy.append(0.)
        if i==5:
            ShapeFuncVectorzz.append(0.)
            ShapeFuncVectorzz.append(0.)
            ShapeFuncVectorzz.append(0.)
            ShapeFuncVectorzz.append(0.)                    
            ShapeFuncVectorzz.append(0.)
            ShapeFuncVectorzz.append(N)                

ShapeFuncVector=np.matrix([ShapeFuncVectorx,ShapeFuncVectory,ShapeFuncVectorz,ShapeFuncVectorxx,ShapeFuncVectoryy,ShapeFuncVectorzz])
print(ShapeFuncVector)
import numpy as np
from sympy import *
import matplotlib as plt

class FEA:
    def __init__(self):
        Check="OK"   
    
    def stiffnessmatrix(self,v=0.3,E=210000):
        Er=(E/((1+v)*(1-(2*v))))
        G=(0.5*E)/(1+v)
        Stiffnessmatrix=np.matrix([[Er*(1.-v),Er*v,Er*v,0.,0.,0.],[Er*v,Er*(1.-v),Er*v,0.,0.,0.],[Er*v,Er*v,Er*(1.-v),0.,0.,0.],[0.,0.,0.,G,0.,0.],[0.,0.,0.,0.,G,0.],[0.,0.,0.,0.,0.,G]])
        #Stiffnessmatrix=np.matrix([[1.,(1./2.)-v,(1./2.)-v],[(1./2.)-v,1.,(1./2.)-v],[(1./2.)-v,(1./2.)-v,1.]])*E

        return Stiffnessmatrix
    
    def shapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        '''
        Shapefunctions are defined for interpolation of displacement inside of the element. To get the strains, the shapefunctions are derived with respect to x, y and z.
        The differentiation has to be done based on a partial basis. first the shapefunctions are differentiated wrp to the natural coordinates (dNexi,...) 
        then the derivates of the actual coordinates wrp to the natural coordinates are carried out. Theese two are then multiplied.
        
        The Jacobian describes how the coordinates of the element relate to the natural coordinates.
        '''

        nodeloc=nodelocations
        dNexi=np.array([-(1-eta)*(1-zeta), (1-eta)*(1-zeta), (1+eta)*(1-zeta),-(1+eta)*(1-zeta),-(1-eta)*(1+zeta), (1-eta)*(1+zeta), (1+eta)*(1+zeta),-(1+eta)*(1+zeta)])*(1/8)
        dNeta=np.array([-(1-exi)*(1-zeta),-(1+exi)*(1-zeta), (1+exi)*(1-zeta), (1-exi)*(1-zeta),-(1-exi)*(1+zeta),-(1+exi)*(1+zeta), (1+exi)*(1+zeta), (1-exi)*(1+zeta)])*(1/8)
        dNzeta=np.array([-(1-exi)*(1-eta), -(1+exi)*(1-eta), -(1+exi)*(1+eta), -(1-exi)*(1+eta),  (1-exi)*(1-eta),  (1+exi)*(1-eta),  (1+exi)*(1+eta),  (1-exi)*(1+eta)])*(1/8)
        J11=J21=J31=J12=J22=J32=J13=J23=J33=0.
      
        for i in range(8):
            n=np.arange(0,25,3)[i]
            xi=nodeloc[n]
            yi=nodeloc[n+1]
            zi=nodeloc[n+2]
            
      
    
            J11=(dNexi[i] *xi)+J11
            J21=(dNeta[i] *xi)+J21
            J31=(dNzeta[i]*xi)+J31
            J12=(dNexi[i] *yi)+J12
            J22=(dNeta[i] *yi)+J22
            J32=(dNzeta[i]*yi)+J32
            J13=(dNexi[i] *zi)+J13
            J23=(dNeta[i] *zi)+J23
            J33=(dNzeta[i]*zi)+J33
        J=np.matrix([[J11,J12,J13],[J21,J22,J23],[J31,J32,J33]])
        Jdet=np.linalg.det(J)
        Jinv=np.linalg.inv(J)

        B=Jinv*np.matrix([dNexi,dNeta,dNzeta])

        Bx=B[0]
        By=B[1]
        Bz=B[2]
  
        B1x= [Bx[0,0],0.,0.,Bx[0,1],0.,0.,Bx[0,2],0.,0.,Bx[0,3],0.,0.,Bx[0,4],0.,0.,Bx[0,5],0.,0.,Bx[0,6],0.,0.,Bx[0,7],0.,0.]
        B2y= [0.,By[0,0],0.,0.,By[0,1],0.,0.,By[0,2],0.,0.,By[0,3],0.,0.,By[0,4],0.,0.,By[0,5],0.,0.,By[0,6],0.,0.,By[0,7],0.]
        B3z= [0.,0.,Bz[0,0],0.,0.,Bz[0,1],0.,0.,Bz[0,2],0.,0.,Bz[0,3],0.,0.,Bz[0,4],0.,0.,Bz[0,5],0.,0.,Bz[0,6],0.,0.,Bz[0,7]]
        B4xy=[By[0,0],Bx[0,0],0.,By[0,1],Bx[0,1],0.,By[0,2],Bx[0,2],0.,By[0,3],Bx[0,3],0.,By[0,4],Bx[0,4],0.,By[0,5],Bx[0,5],0.,By[0,6],Bx[0,6],0.,By[0,7],Bx[0,7],0.]
        B5zy=[0.,Bz[0,0],By[0,0],0.,Bz[0,1],By[0,1],0.,Bz[0,2],By[0,2],0.,Bz[0,3],By[0,3],0.,Bz[0,4],By[0,4],0.,Bz[0,5],By[0,5],0.,Bz[0,6],By[0,6],0.,Bz[0,7],By[0,7]]
        B6zx=[Bz[0,0],0.,Bx[0,0],Bz[0,1],0.,Bx[0,1],Bz[0,2],0.,Bx[0,2],Bz[0,3],0.,Bx[0,3],Bz[0,4],0.,Bx[0,4],Bz[0,5],0.,Bx[0,5],Bz[0,6],0.,Bx[0,6],Bz[0,7],0.,Bx[0,7]]
        Btot=np.matrix([B1x,B2y,B3z,B4xy,B5zy,B6zx])
           
        return Btot, Jdet
             
    def Volume(self,nodelocations=[]):
        
        n=nodelocations
        tetrahedrons=[[1,6,8,5],[1,4,6,2],[1,6,8,4],[7,6,8,3],[3,4,6,2],[8,7,6,4]]
        Volume=0.
        for tet in tetrahedrons:
        
            t=np.array(tet)-[1,1,1,1]
            a1=n[t[0]*3]
            a2=n[(t[0]*3)+1]
            a3=n[(t[0]*3)+2]
            b1=n[(t[1]*3)]
            b2=n[(t[1]*3)+1]
            b3=n[(t[1]*3)+2]
            c1=n[(t[2]*3)]
            c2=n[(t[2]*3)+1]
            c3=n[(t[2]*3)+2]
            d1=n[(t[3]*3)]
            d2=n[(t[3]*3)+1]
            d3=n[(t[3]*3)+2]
            Volume=((1/6.)*abs(np.linalg.det(np.matrix([[a1-d1,b1-d1,c1-d1],[a2-d2,b2-d2,c2-d2],[a3-d3,b3-d3,c3-d3]]))))+Volume    

        return Volume
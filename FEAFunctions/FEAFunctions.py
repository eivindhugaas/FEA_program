import numpy as np
from sympy import *
import matplotlib as plt

class FEA:
    def __init__(self):
        Check="OK"    
    def shapefunction8node(self,xi=0.5,eta=1.0,zeta=0.1):
        '''
        isoparametric, hexahedral, natural coordinates:
        xi is an e with a krussedull on top
        eta is an n with a long leg
        zeta is a snake that is bowing, it is zzzzzzeeeetaaaaaa.
        '''
        
        nodecoords=[[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.],[1.,-1.,1.],[1.,1.,1.],[-1.,1.,1.]]
        ShapeFuncVectorx=[]
        ShapeFuncVectory=[]
        ShapeFuncVectorz=[]
        ShapeFuncVector=[]
        for i in range(0,3,1):
            for coord in nodecoords:
                x=coord[0]
                y=coord[1]
                z=coord[2]
                N=(1/8.)*(1.+(xi*x))*(1.+(eta*y))*(1.+(zeta*z))*((xi*x)+(eta*y)+(zeta*z)-2.)
                if i==0:
                    ShapeFuncVectorx.append(N)
                    ShapeFuncVectorx.append(0.)
                    ShapeFuncVectorx.append(0.)
                if i==1:
                    ShapeFuncVectory.append(0.)                    
                    ShapeFuncVectory.append(N)
                    ShapeFuncVectory.append(0.)
                if i==2:
                    ShapeFuncVectorz.append(0.)                    
                    ShapeFuncVectorz.append(0.)
                    ShapeFuncVectorz.append(N)

        ShapeFuncVector.append([ShapeFuncVectorx,ShapeFuncVectory,ShapeFuncVectorz])
        ShapeFuncVector=np.matrix(ShapeFuncVector[0])
        return ShapeFuncVector
    
    #def element8node(self,disp=[],location=[],xi=0.5,eta=1.0,zeta=0.1):
        #'''
        #Disp must be in format=[..........] and is an array with three disp values for each node in successive order.
        #Location must be in same format as disp and is the physical location of the nodes.
        #'''
        #newlocations=np.transpose(np.add(np.matrix(disp),np.matrix(location)))
        #shapefunc=FEA.shapefunction8node(self,xi=xi,eta=eta,zeta=zeta)
        #locationinpoint=np.dot(shapefunc,newlocations)
        #ex=(1/4.)*(disp[]
        #elementstrain=[ex,ey,ez,txy,txz,tyz]
        
        #return locationinpoint
    
    def stiffnessmatrix(self,v=0.3,E=210000):
        Er=(E/((1+v)*(1-(2*v))))
        G=(0.5*E)/(1+v)
        Stiffnessmatrix=np.matrix([[Er*(1.-v),Er*v,Er*v,0.,0.,0.],[Er*v,Er*(1.-v),Er*v,0.,0.,0.],[Er*v,Er*v,Er*(1.-v),0.,0.,0.],[0.,0.,0.,G,0.,0.],[0.,0.,0.,0.,G,0.],[0.,0.,0.,0.,0.,G]])
        #Stiffnessmatrix=np.matrix([[1.,(1./2.)-v,(1./2.)-v],[(1./2.)-v,1.,(1./2.)-v],[(1./2.)-v,(1./2.)-v,1.]])*E

        return Stiffnessmatrix
  
    def HydrostatStresshapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        N=1.        
        return N
    
    def AvgStresshapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        N=(1./8.)*np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])      
        Amat=[]
        ver=[]
        for t in range(0,6):
            hor=[]
            for i in range(0,48):
                a=0.
                if i==0+t or i==6+t or i==12+t or i==18+t or i==24+t or i==30+t or i==36+t or i==42+t:
                    a=1./8.
                hor.append(a)
            Amat.append(hor)
        N=np.matrix(Amat)  
        N=np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])  
        return N
    
    def AvgStrainshapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        N=(1./8.)*np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
        Amat=[]
        ver=[]
        for t in range(0,6):
            hor=[]
            for i in range(0,48):
                a=0.
                if i==0+t or i==6+t or i==12+t or i==18+t or i==24+t or i==30+t or i==36+t or i==42+t:
                    a=1./8.
                hor.append(a)
            Amat.append(hor)
        N=np.matrix(Amat)  
        


        N=np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])  
        return N   
    
    def UStraintoAvgStrainshapefunc(self):

        Amat=[]
        ver=[]
        for t in range(0,6):
            hor=[]
            for i in range(0,6):
                a=0.
                if i==0+t or i==6+t or i==12+t or i==18+t or i==24+t or i==30+t or i==36+t or i==42+t:
                    a=1./8.
                hor.append(a)
            Amat.append(hor)
        N=np.matrix(Amat)
       
        return N     
    
    def norderStrainshapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        '''
        isoparametric, hexahedral, natural coordinates:
        xi is an e with a krussedull on top
        eta is an n with a long leg
        zeta is a snake that is bowing, it is zzzzzzeeeetaaaaaa.
        '''
        
        nodecoords=[[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.],[1.,-1.,1.],[1.,1.,1.],[-1.,1.,1.]] #in natural coords
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
                #N=(1/8.)*(1.+(exi*x))*(1.+(eta*y))*(1.+(zeta*z))*((exi*x)+(eta*y)+(zeta*z)-2.)
                N=(1./8.)*(1.+(x*exi))*(1.+(eta*y))*(1+(zeta*z))

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
        
        return ShapeFuncVector         
    
    def twoorderStrainshapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        '''
        isoparametric, hexahedral, natural coordinates:
        xi is an e with a krussedull on top
        eta is an n with a long leg
        zeta is a snake that is bowing, it is zzzzzzeeeetaaaaaa.
        '''
        
        nodecoords=[[-1.,-1.,-1.],[1.,-1.,-1.],[1.,1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.],[1.,-1.,1.],[1.,1.,1.],[-1.,1.,1.]] #in natural coords
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
                #N=(1/8.)*(1.+(exi*x))*(1.+(eta*y))*(1.+(zeta*z))*((exi*x)+(eta*y)+(zeta*z)-2.)
                N=(1./8.)*(1.+(x*exi))*(1.+(eta*y))*(1+(zeta*z))

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
        
        return ShapeFuncVector             

    
    def shapefunc(self,nodelocations=[],exi=0.5,eta=1.0,zeta=0.1):
        '''
        Shapefunctions are defined for interpolation of displacement inside of the element. To get the strains, the shapefunctions are derived with respect to x, y and z.
        The differentiation has to be done based on a partial basis. first the shapefunctions are differentiated wrp to the natural coordinates (dNexi,...) 
        then the derivates of the actual coordinates wrp to the natural coordinates are carried out. Theese two are the multiplied.
        
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

        B=np.matmul(Jinv,((np.matrix([dNexi,dNeta,dNzeta]))/Jdet))
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
    
    
    
    def plotelement(self,nodelocations=[]):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        r = [-1,1]
        
        X, Y = np.meshgrid(r, r)
        
        ax.plot_surface(X,Y,1, alpha=0.5)
        ax.plot_surface(X,Y,-1, alpha=0.5)
        ax.plot_surface(X,-1,Y, alpha=0.5)
        ax.plot_surface(X,1,Y, alpha=0.5)
        ax.plot_surface(1,X,Y, alpha=0.5)
        ax.plot_surface(-1,X,Y, alpha=0.5)
        
        ax.scatter3D(points[:, 0], points[:, 1], points[:, 2])
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        plt.show()        
        


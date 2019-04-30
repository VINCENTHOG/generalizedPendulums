from sympy import *
import numpy as np
import time
global g
g = symbols("g")

class Lagrange():
    def __init__(self,nodes):
        self.nodes = nodes
        self.timer = time.time()
    
    def updateTimer(self):
        t          = np.round(time.time()-self.timer,3)
        self.timer = time.time()
        return t

    # Lagrangian's Jacobian with respect to q and dq
    def derivFirst(self,L,q,dq):
        ldq = []
        lq  = []
        for _q,_dq in zip(q,dq):
            lq.append(diff(L,_q))   # Jacobian with respect to vector q
            ldq.append(diff(L,_dq)) # Jacobian with respect to vector dq
        print("\n-- First Derivation:      ",self.updateTimer()," seconds")
        return lq,ldq

    # Time derivative of the Lagrangian's Jacobian with respect to dq
    def derivSecond(self,ldq,q,dq,ddq):
        ldqdt = []
        for _ldq in ldq:
            deriv=0
            for _q,_dq,_ddq in zip(q,dq,ddq):
                deriv+=(diff(_ldq,_q)*_dq+diff(_ldq,_dq)*_ddq) #Time derivative 
            ldqdt.append(deriv)
        print("-- Second Derivation:     ",self.updateTimer()," seconds")
        return ldqdt

    # Substract the two left term in Euler-Lagrange equation
    def computeEulerLagrange(self,ldqdt,lq):
        links = []
        for _ldqdt,_lq in zip(ldqdt,lq):
            links.append(_ldqdt-_lq)
        print("-- Links Creation:        ",self.updateTimer()," seconds")
        return links
    
    # This function is a little trick to isolate variables in Euler-Lagrange formula
    def createCoefficients(self,links,q,dq,ddq):
        global g
        # Generate a polynomial with Sympy
        pddq,pdq,pg = [],[],[] #Polynomial
        cddq,cdq,cg = [],[],[] #Coefficients
        for link in links:
            pddq.append(Poly(link,ddq))
            pdq.append(Poly(link,dq))
            pg.append(Poly(link,g))
        # Isolate polynomial's coefficients
        for _pddq,_pdq,_pg in zip(pddq,pdq,pg):
            cddq.append(_pddq.coeffs())
            cdq.append(_pdq.coeffs())
            cg.append(_pg.coeffs())
        print("-- Coefficients Creation: ",self.updateTimer()," seconds")
        return cddq,cdq,cg
    
    # Generalized Inertia Matrix
    def createMatrix_D(self,cddq):
        D = zeros(len(cddq))
        for i in range(len(cddq)):
            childs = self.nodes[i].childs
            nb_anc = len(self.nodes[i].ancest)
            offset = nb_anc
            for j in range(i+1):
                if(j<i):
                    D[i,j] = D[j,i]
                elif(j==i):
                    D[i,j] = cddq[i][offset]
                    offset += 1
            for k in range(len(childs)):
                D[i,childs[k]] = cddq[i][offset+k]
        print("-- D Matrix Creation:     ",self.updateTimer()," seconds")
        return D

    # Centrifugal and Coriolis Matrix  
    def createMatrix_H(self,cdq):
        H = zeros(len(cdq))
        for i in range(len(cdq)):
            childs = self.nodes[i].childs
            nb_anc = len(self.nodes[i].ancest)
            offset = nb_anc
            for j in range(i+1):
                if(j<i):
                    H[i,j] = -H[j,i]
                elif(j==i):
                    H[i,j] = 0
            for k in range(len(childs)):
                H[i,childs[k]] = cdq[i][offset+k]
        print("-- H Matrix Creation:     ",self.updateTimer()," seconds")
        return H
    
    # Centrifugal and Coriolis Matrix * Speed
    def multiplyDQ(self,H,dq):
        square_dq = Matrix(np.asarray(dq)**2)
        hdq       = H*square_dq
        return hdq
    
    # Generalized Gravitational Matrix
    def createMatrix_G(self,cg):
        global g
        G = zeros(len(cg),1)
        for i in range(len(cg)):
            G[i] = cg[i][0]*g
        print("-- G Matrix Creation:     ",self.updateTimer()," seconds")
        return G
    
    # General extended form of Euler-Lagrange Formula
    #  
    #   D(q)*ddq + H(q,dq) + G(q) = J^t*External_Forces + Torques
    #
    # Right hand side can be completely ignored if not needed 
    def euler_lagrange(self,k,p,q,dq,ddq):
        L           = k-p
        lq,ldq      = self.derivFirst(L,q,dq)                   # Lagrangian's Jacobian with respect to q and dq       
        ldqdt       = self.derivSecond(ldq,q,dq,ddq)            # Time derivation of the Lagrangian's Jacobian with respect to dq
        links       = self.computeEulerLagrange(ldqdt,lq)       # Computation of the difference in the Euler-Lagrange formula
        cddq,cdq,cg = self.createCoefficients(links,q,dq,ddq)   # Trick to isolate q, dq and ddq
        D           = self.createMatrix_D(cddq)                 # Generalized Inertia Matrix
        H           = self.createMatrix_H(cdq)                  # Centrifugal and Coriolis Matrix
        Hdq         = self.multiplyDQ(H,dq)                     # Computation of the dot product between Centrifugal and Coriolis Matrix and dq
        G           = self.createMatrix_G(cg)                   # Generalized Gravitational Matrix
        return D,Hdq,G
    
'''
    EXTRA DETAILS
   
    Kinetic Energy:

            K = 1/2 * mass * speed**2
            ** Speed != dq

    Potential Energy:

            P = mass * gravitational acceleration * height

    Lagrangian:
    
            Lagrangian = Kinetic Energy - Potential Energy
                L      =    K           -       P    
    
    Euler-Lagrange Formula:
    
            d(del(L)/del(dq))/dt - del(L)/del(q) = 0
    
    in other words...
    
            time_derivative(Lagrangian's Jacobian with respect to dq) - Lagrangian's Jacobian with respect to q = 0
    
    Isolate ddq and dq
    
            D(q)*ddq 
            and
            H(q)*dq
    
    The rest becomes...
    
            G(q)
    
    Final Result
    
            D(q)*ddq + H(q,dq) + G(q) = 0
'''

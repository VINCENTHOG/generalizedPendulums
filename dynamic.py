
import numpy as np
from numpy import sin,cos
import sys,os,parser
DIMENSIONS = 2
RESIZE_STR = 7

class Dynamic:
    def __init__(self,D,H,G,J,K,P):
        self.D = parser.expr(str(D)[RESIZE_STR:-1]).compile()
        self.H = parser.expr(str(H)[RESIZE_STR:-1]).compile()
        self.G = parser.expr(str(G)[RESIZE_STR:-1]).compile()
        self.J = parser.expr(str(J)[RESIZE_STR:-1]).compile()
        self.K = parser.expr(str(K)).compile()
        self.P = parser.expr(str(P)).compile()
        
    def updateMatrix(self,param):
        g = 9.81
        q = param[0:int(len(param)/2)];dq=param[int(len(param)/2):len(param)] 
        d = eval(self.D)
        h = eval(self.H)
        g = eval(self.G)
        j = eval(self.J)
        return np.asmatrix(d),np.asmatrix(h),np.asmatrix(g),np.asmatrix(j)

    # Solving acceleration in equation of motion
    # ddq = inverse(D(q)) * (F + U - H(q,dq) - G(q))
    def deriv(self,param,u,fext):
        U       = np.array([u]).T
        D,H,G,J = self.updateMatrix(param)
        fext    = np.array([fext]).T
        F       = J * fext
        fuhg    = F + U - H - G
        inv_D   = D.I.A
        dq      = np.array(param[int(len(param)/2):len(param)])
        ddq     = np.dot(inv_D,fuhg)
        return np.append(dq,ddq)

    # One step in environment
    def step(self,param,timeStep,u,fext):
        # Little Runge-Kutta 4th order for better derivative estimation
        _k1    = self.deriv(param,u,fext)
        _k2    = self.deriv(param + 0.5 * timeStep * _k1,u,fext)
        _k3    = self.deriv(param + 0.5 * timeStep * _k2,u,fext)
        _k4    = self.deriv(param + timeStep*_k3,u,fext)
        param  = param + (timeStep/6.0) * (_k1 + 2 * _k2 + 2 * _k3 + _k4)
        energy = self.getMechanicalEnergy(param)
        print("System Energy = ",energy)
        return param, energy

    # Kinetic and Potential Energy
    def getMechanicalEnergy(self,param):
        try:
            g  = 9.81
            q  = param[0:int(len(param)/2)]
            dq = param[int(len(param)/2):len(param)]
            k  = eval(self.K)
            p  = eval(self.P)
            return k + p
        except:
            print("\n\t--Can't display Energy: Wrong Model\n\t--Change getMechanicalEnergy() equations\n")
            return 0


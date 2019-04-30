
import numpy as np
from numpy import sin,cos
import sys,os
DIMENSIONS = 2

# D,H,G and J need to be computed before simulation
# $ python model.py 1
# Copy paste them under
def updateMatrix(param):
    try:
        g = 9.81
        q = param[0:int(len(param)/2)];dq=param[int(len(param)/2):len(param)] 
        D = np.matrix([[1.58333333333333, 0.5*cos(q[0] - q[1])], [0.5*cos(q[0] - q[1]), 0.583333333333333]])
        H = np.matrix([[0.5*dq[1]**2*sin(q[0] - q[1])], [-0.5*dq[0]**2*sin(q[0] - q[1])]])
        G = np.matrix([[1.5*g*cos(q[0])], [0.5*g*cos(q[1])]])
        J = np.matrix([[-1.0*sin(q[0]), 1.0*cos(q[0]), -1.0*sin(q[0]), 1.0*cos(q[0])], [0, 0, -1.0*sin(q[1]), 1.0*cos(q[1])]])
        return D,H,G,J
    except:
        print("\n--Error in the equation of motion. Make sure the right model is loaded\n\nExit Program.")
        os._exit(0)

# Solving acceleration in equation of motion
# ddq = inverse(D(q)) * (F + U - H(q,dq) - G(q))
def deriv(param,u,fext):
    U       = np.array([u]).T
    D,H,G,J = updateMatrix(param)
    fext    = np.array([fext]).T
    F       = J * fext
    fuhg    = F + U - H - G
    inv_D   = D.I.A
    dq      = np.array(param[int(len(param)/2):len(param)])
    ddq     = np.dot(inv_D,fuhg)
    return np.append(dq,ddq)

# One step in environment
def step(param,timeStep,u,fext):
    # Little Runge-Kutta 4th order for better derivative estimation
    _k1    = deriv(param,u,fext)
    _k2    = deriv(param + 0.5 * timeStep * _k1,u,fext)
    _k3    = deriv(param + 0.5 * timeStep * _k2,u,fext)
    _k4    = deriv(param + timeStep*_k3,u,fext)
    param  = param + (timeStep/6.0) * (_k1 + 2 * _k2 + 2 * _k3 + _k4)
    energy = getMechanicalEnergy(param)
    print("System Energy = ",energy)
    return param, energy

# Kinetic and Potential Energy
def getMechanicalEnergy(param):
    # Mechanical energy appears when you run:
    #       python model.py 1
    # Copy paste the result under
    try:
        g  = 9.81
        q  = param[0:int(len(param)/2)]
        dq = param[int(len(param)/2):len(param)]
        K = 0.125*dq[0]**2*sin(q[0])**2 + 0.125*dq[0]**2*cos(q[0])**2 + 0.166666666666667*dq[0]**2 + 0.166666666666667*dq[1]**2 + 0.5*(-1.0*dq[0]*sin(q[0]) - 0.5*dq[1]*sin(q[1]))**2 + 0.5*(1.0*dq[0]*cos(q[0]) + 0.5*dq[1]*cos(q[1]))**2
        P = 1.0*g*(1.0*sin(q[0]) + 0.5*sin(q[1])) + 0.5*g*sin(q[0])
        return K + P
    except:
        print("\n\t--Can't display Energy: Wrong Model\n\t--Change getMechanicalEnergy() equations\n")
        return 0


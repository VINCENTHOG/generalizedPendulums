from sympy import *
from node import Node
from lagrange import Lagrange
from pendulums import *
import numpy as np
import sys,os

DEFAULT_LENGTH    = 1.0
DEFAULT_MASS      = 1.0
DEFAULT_COM_RATIO = 0.5
MAX_NODE          = 40
MOTHER_ID         = 777         #Don't touch :)
EQUATION          = 1
DISPLAY           = 0

global  g
g  = symbols("g")

class Model():
    def __init__(self):
        self.nodes        = []
        self.node_counter = 0
        self.q            = []
        self.dq           = []
        self.ddq          = []

    # Creates a new node in the tree with a specific parent
    def createNode(self,l,m,parent,com):
        if(self.node_counter<MAX_NODE):
            q       = symbols("q["+str(self.node_counter)+"]")
            inertia = 1/3*m*l**2
            n       = Node(self.node_counter,l,m,inertia,q,com)
            n.setPos(parent)
            n.set_ancest()
            self.nodes.append(n)
            self.node_counter+=1
        else:
            print("Maximum node reached")

    # Sympy symbols creation
    def createVariables(self):
        for node in self.nodes:
            self.q.append(symbols("q["+str(node.id)+"]"))
            self.dq.append(symbols("dq["+str(node.id)+"]"))
            self.ddq.append(symbols("ddq["+str(node.id)+"]"))

    # Time derivative of the global positions 
    def derivModel(self):
        global nodes
        for i in range(self.node_counter):
            _dx=0;_dy=0
            for j in range(i+1):
                _dx+=diff(self.nodes[i].com_x,self.q[j])*self.dq[j]
                _dy+=diff(self.nodes[i].com_y,self.q[j])*self.dq[j]
            self.nodes[i].dx = _dx
            self.nodes[i].dy = _dy
    
    # Computes Potential and Kinetic energies
    def getEnergies(self):
        global g
        K,P = 0,0
        for i in range(self.node_counter):
            K+=0.5*self.nodes[i].m*(self.nodes[i].dx**2+self.nodes[i].dy**2)+0.5*self.nodes[i].i*self.dq[i]**2
            P+=self.nodes[i].m*g*self.nodes[i].com_y
        print("Kinetic and Potential Energy:\n\nK =",K,"\nP =",P)
        return K,P

    # Little sorting method
    def sortAncesters(self,ancesters):
        count = 0
        for i in range(len(ancesters)-1):
            count+=1
            tmp   = ancesters[i+1]
            if(len(ancesters[i])<len(tmp)):
                j=1
                ancesters[i+1] = ancesters[i]
                ancesters[i]   = tmp
                while(not(i-j<0) and len(ancesters[i-j])<len(tmp)):
                    ancesters[i-j+1] = ancesters[i-j]
                    ancesters[i-j]   = tmp
                    j+=1
                    count+=1
        return ancesters

    # Defines tree linkage while avoiding going back to root every time
    def createLines(self,ancesters):
        lines = []
        for ancester in ancesters:
            if(len(lines)==0):
                lines.append(ancester)
            else:
                index = len(ancester) - 1
                new   = True
                for line in lines:
                    if (line[index]==ancester[index]):
                        new = False
                        break
                if(new):
                    lines.append(ancester)
        return lines

    # Call to the sympy simplify method. 
    # CAUTION: AVOID TO CALL THIS FUNCTION IF THERE IS MORE THAN 10 NODES***
    def simplifyLagrangian(self,D,H,G):
        print("\nSimplification...\n\n\n-- Matrix D:\n",simplify(D),"\n\nMatrix H:\n",simplify(H),"\n\nMatrix G:\n",simplify(G))
    
    # Prints the jacobian of the generalized position matrix
    def jacobianGeneralizedPositionMatrix(self):
        model = Matrix([[]])
        for node in self.nodes:
            model = model.row_insert(model.shape[0],Matrix([[node.pos_x]]))
            model = model.row_insert(model.shape[0],Matrix([[node.pos_y]]))
        print("\nTranspose Jacobian of phi(q) ")  
        print(simplify(model.jacobian(self.q).T))

    # Creates your model
    # Use this method to add/remove nodes
    def generateModel(self):
        l,m,i                = 0,0,0
        symbol               = symbols("mother")
        mother_of_all        = Node(MOTHER_ID,l,m,i,symbol,DEFAULT_COM_RATIO)
        mother_of_all.parent = mother_of_all
        parent               = mother_of_all
        try:
            for i in range(2):
               self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,parent,DEFAULT_COM_RATIO)
               parent = self.nodes[-1]
            # FOLLOWING LINES ARE AN EXAMPLE. CREATE YOUR OWN :)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[0],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[3],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[5],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[0],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[10],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[13],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[15],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[0],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[8],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[17],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[19],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[6],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[7],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[10],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[21],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[14],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[28],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[25],DEFAULT_COM_RATIO)
            # self.createNode(DEFAULT_LENGTH,DEFAULT_MASS,self.nodes[-1],DEFAULT_COM_RATIO)
        except:
            print("\nError in your nodes\nExit Program.")
            os._exit(0)
    
    #Calls the Euler-Lagrange generator from lagrange.py
    def generateEulerLagrange(self):
        self.createVariables()
        self.derivModel()
        k,p   = self.getEnergies()
        L     = Lagrange(self.nodes)
        D,H,G = L.euler_lagrange(k,p,self.q,self.dq,self.ddq) # H is actually H*(dq**2)
        try: 
            answer = input("Do you want to simplify the lagrangian?(Might take some time) \n\n[y/n]:")
            if(answer == "y" or answer =="Y"):
                self.simplifyLagrangian(D,H,G)
            else:
                print("\n\n-- Matrix D:\n",D,"\n\nMatrix H:\n",H,"\n\nMatrix G:\n",G)
        except:
            print("Don't know what you did but...")
            print("\n\n-- Matrix D:\n",D,"\n\nMatrix H:\n",H,"\n\nMatrix G:\n",G)
        self.jacobianGeneralizedPositionMatrix()
    
    # Matplotlib display of your model 
    # CAUTION: Make sure that dynamic.py has the proper dynamic and mechanical energy 
    def displayModel(self):
        ancest_table = []
        for node in self.nodes:
            ancest_table.append(np.hstack((node.ancest,int(node.id))))
        ancest_table = self.sortAncesters(ancest_table)
        lines        = self.createLines(ancest_table)
        lengths      = [n.l for n in self.nodes]
        pendulums    = Pendulums(self.node_counter,lines,lengths)
        pendulums.run()
        pendulums.displayEnergy()
    
    # Argument parser
    def runScenario(self):
        try:
            scenario = int(sys.argv[1])
            if (scenario == DISPLAY or scenario == EQUATION):
                if(scenario == DISPLAY):
                    self.displayModel()
                elif(scenario == EQUATION):
                    self.generateEulerLagrange()
            else:
                print("Invalid Argument\n\n0: Simulation\n1: Equations Computation")
        except:
            print("\n--No Argument--\n\n0: Simulation\n1: Equations Computation")
        

if __name__ == "__main__":
    m = Model()
    m.generateModel()
    m.runScenario()

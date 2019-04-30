import numpy as np
from sympy import *

MOTHER_ID         = 777 #Don't touch :)

class Node():
    
    def __init__(self,_id,l,m,i,symbol,com):
        self.parent = 0                                 #Parent Node
        self.l      = l                                 #Length
        self.m      = m                                 #Mass
        self.i      = i                                 #Inertia
        self.com_r  = com                               #Ratio of the distance of center of mass from connection point
        self.id     = _id                               #Node ID(Chronological)
        self.q      = symbol                            #Sympy symbol
        self.pos_x  = self.l*cos(self.q)                #Position of the termination of the pendulum (X)
        self.pos_y  = self.l*sin(self.q)                #Position of the termination of the pendulum (Y)
        self.com_x  = self.l*self.com_r*cos(self.q)     #Position of the center of mass of the pendulum (X)
        self.com_y  = self.l*self.com_r*sin(self.q)     #Position of the center of mass of the pendulum (X)
        self.ancest = []                                #All the previous nodes till root
        self.childs = []                                #All the next node till leaves
        self.dx     = 0                                 #Velocity (X)
        self.dy     = 0                                 #Velocity (Y)

    #Setting global position of center of mass, termination and the parent node
    def setPos(self,node):
        self.parent = node      
        self.pos_x  += self.parent.pos_x 
        self.pos_y  += self.parent.pos_y 
        self.com_x  += self.parent.pos_x 
        self.com_y  += self.parent.pos_y
    
    #Setting the ancesters list until root from a given parent node
    def set_ancest(self):
        if(not(self.parent==0)):
            p = self.parent
            while(not(p.id==MOTHER_ID)):
                p.childs.append(self.id)
                self.ancest.append(p.id)
                p = p.parent
        self.ancest.reverse()

    def show(self):
        print("\n------ id:",self.id,"\np: ",self.parent.id,"\nl: ",self.l,"\nm: ",self.m,"\ni: ",self.i,"\npx:",self.pos_x,"\npy:",self.pos_y,"\ncx:",self.com_x,"\ncy:",self.com_y,"\nan:",self.ancest,"\nch:",self.childs)


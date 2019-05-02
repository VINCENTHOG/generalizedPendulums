import numpy as np
import dynamic as dynamic
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from numpy import array
import sys

DIMENSIONS = 2

class Pendulums():
    def __init__(self,nb,lines,lengths):
        self.nb       = nb
        self.q        = np.random.rand(self.nb)*np.pi
        self.dq       = np.zeros(self.nb)
        self.param    = np.hstack((self.q,self.dq))
        self.lines    = lines
        self.lengths  = lengths
        self.min_x    = -len(lines[0]) - 0.5            #Viewer
        self.min_y    = -len(lines[0]) - 0.5            #Viewer
        self.max_x    =  len(lines[0]) + 0.5            #Viewer
        self.max_y    =  len(lines[0]) + 0.5            #Viewer
        self.fext     = np.zeros(nb*DIMENSIONS)
        self.torque   = np.zeros(self.nb)
        self.timeStep = 0.01
        self.energy   = []

    def updateParam(self,param,energy):
        self.param = param
        self.q     = param[0:self.nb]
        self.dq    = param[self.nb:len(self.param)]
        self.fext  = np.zeros(self.nb*DIMENSIONS)
        self.energy.append(energy)

    def displayEnergy(self):
        fig = plt.figure()
        fig.suptitle('Mechanical Energy', fontsize=16)
        ax  = fig.add_subplot(111,frameon=True, autoscale_on=True)
        ax.plot(np.arange(len(self.energy)),self.energy,  lw=2,label="Energy")
        ax.set_xlabel('Iterations')
        ax.set_ylabel('Mechanical Energy')
        plt.show()

    def updateLinks(self):
        aline_x = []
        aline_y = []    
        for line in self.lines:
            mline_x = [0.0]
            mline_y = [0.0]
            for i in line:
                mline_x.append(mline_x[-1] + self.lengths[int(i)] * np.cos(self.q[int(i)]))
                mline_y.append(mline_y[-1] + self.lengths[int(i)] * np.sin(self.q[int(i)]))
            aline_x.append(mline_x)
            aline_y.append(mline_y)
        return aline_x,aline_y

    def run(self):
        fig    = plt.figure(figsize=(10, 10))
        ax     = fig.add_subplot(111,frameon=True, autoscale_on=False, xlim=(self.min_x, self.max_x), ylim=(self.min_y, self.max_y))
        lines  = []
        for i in range(self.nb):
            line = ax.plot([], [], '-k', lw=2)[0]
            lines.append(line)

        def init():
            for line in lines:
                line.set_data([], [])            
            return lines

        def animate(i):
            param,energy = dynamic.step(self.param,self.timeStep,self.torque,self.fext)
            self.updateParam(param,energy)
            x,y = self.updateLinks()
            for _x,_y,line in zip(x,y,lines):
                line.set_data(_x,_y)
            return lines

        ani = animation.FuncAnimation(fig, animate, np.arange(0,100000),interval=10, blit=True, init_func=init,repeat_delay=100)
        # Next line saves the animation (Slows down computation)
        # ani.save('pendulum.mp4', fps=100)
        plt.show()
    

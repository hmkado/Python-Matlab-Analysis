# -*- coding: utf-8 -*-
"""

@author: JR
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.physics.mechanics import *

#4.1
print('#4.1')
m,J,k1,c1,k2,c2,l1,l2 = sp.symbols('m J k1 c1 k2 c2 l1 l2', real= True)
x = dynamicsymbols('x')
xdot = dynamicsymbols('x', 1)
theta = dynamicsymbols('theta')
thetadot = dynamicsymbols('theta', 1)

T = (1/2)*J*(thetadot)**2+(1/2)*m*(xdot)**2
U = (1/2)*k1*(x+l1*theta)**2+(1/2)*k2*(x-l2*theta)**2
L = T-U

#Eqn 1
eq1 = sp.simplify(sp.diff(sp.diff(L, xdot), 't') - sp.diff(L, x))

#Eqn 2
eq2 = sp.simplify(sp.diff(sp.diff(L, thetadot), 't') - sp.diff(L, theta))

#damping
Fe1 = -c1*(xdot+l1*thetadot)
Fe2 = -c2*(xdot-l2*thetadot)
print('EoM:')
print('In terms of x: ',eq1,'=',Fe1)
print('In terms of theta: ',eq2,'=',Fe2)

#4.2
print('\n#4.2')
M = np.array([[m, 0],[0, J]])
C = np.array([[-c1, -c1*l1],[-c2, c2*l2]])
K = np.array([[k1+k2, k1*l1-k2*l2],[k1-k2, k1*l1+k2*l2]])

xdot2 = np.transpose(np.array([sp.diff(xdot), sp.diff(thetadot)]))
xdotm = np.transpose(np.array([xdot, thetadot]))
xm = np.transpose(np.array([x, theta]))
EoMm = np.matmul(M,xdot2)+np.matmul(C,xdotm)+np.matmul(K,xm)
print('EoM: ',EoMm)

#4.3
print('\n#4.3')
x0,xdot0=(.1,-1)
theta0,thetadot0=(.05,.3)
m,J=(1200,300)
k1=k2=10*10**3
c1=c2=0
l1,l2=(1.5,.8)
print(l2)

n,tf=(100000,10)
t=np.zeros(n+1)
t[0]=0
delt=(tf-t[0])/n

xnum=np.zeros(n+1)
dxnum=np.zeros(n+1)
xnum[0]=x0
dxnum[0]=xdot0

thetanum=np.zeros(n+1)
dthetanum=np.zeros(n+1)
thetanum[0]=theta0
dthetanum[0]=thetadot0

for i in range(n):
    
    a1=-(k1*(xnum[i]+l1*thetanum[i])+k2*(xnum[i]-l2*thetanum[i]))/m
    
    xnum[i+1]=xnum[i]+dxnum[i]*delt
    dxnum[i+1]=dxnum[i]+a1*delt
    
    a2=-(k1*(xnum[i]+l1*thetanum[i])-k2*(xnum[i]-l2*thetanum[i]))/J
    
    thetanum[i+1]=thetanum[i]+dthetanum[i]*delt
    dthetanum[i+1]=dthetanum[i]+a2*delt
    
    t[i+1]=t[i]+delt

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)

ax1.plot(t, xnum)
ax1.set_title('x [m]')
ax2.plot(t, thetanum)
ax2.set_title('Theta [rad]')
ax3.plot(t, dxnum)
ax3.set_title('x dot [m/s]')
ax4.plot(t, dthetanum)
ax4.set_title('Theta dot [rad/s]')

plt.savefig('fig.png')
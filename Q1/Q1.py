# -*- coding: utf-8 -*-
"""

@author: JR
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.physics.mechanics import *

#1.1
print('#1.1')
m,r,k,c,J = sp.symbols('m r k c J', real= True)
x = dynamicsymbols('x')
xdot = dynamicsymbols('x', 1)
T=(1/2)*(J*(xdot/r)**2)+(1/2)*m*(xdot**2)
U=(1/2)*k*x**2
L=T-U
LHS = sp.simplify(sp.diff(sp.diff(L, xdot), 't') - sp.diff(L, x)) #d/dt(dL/dxdot)-dL/dx
RHS = -c*xdot
f=LHS-RHS
print('Equation of motion = ',f,' = ',0)

#1.2
print('\n#1.2')
k=10*(10**3)
m=12
r=.34
J=.016
c=2*np.sqrt(k*m)
print('dashpot damping ratio less than ',c,'N/m/s keeps damping ratio less than 1')

#1.3
print('\n#1.3')
c=31
x0=.03
v0=-.2
n=100000
t=np.zeros(n+2)
tf=5
t[1]=0
delt=(tf-t[1])/n

x=np.zeros(n+2)
v=np.zeros(n+2)
x[1]=x0
v[1]=v0
for i in range(n):
    x[i+2]=x[i+1]+v[i+1]*delt
    v[i+2]=v[i+1]+((-c/((J/r**2)+m))*v[i+1] + (-k/((J/r**2)+m))*x[i+1])*delt
    t[i+2]=t[i+1]+delt
    
fig, (ax1, ax2) = plt.subplots(2)

ax1.plot(t, x)
ax1.set_title('Position [m]')
ax1.set_xlim([0, 5])
ax1.set_ylim([-.04, .04])
ax1.axhline(y=.01,linestyle='dashed')
ax1.axhline(y=-.01,linestyle='dashed')
ax1.axvline(x=1,linestyle='dashed')
ax2.plot(t, v)
ax2.set_title('Velocity [m/s]')
ax2.set_xlim([0, 5])
ax2.set_ylim([-1, 1])
ax2.axhline(y=.3,linestyle='dashed')
ax2.axhline(y=-.3,linestyle='dashed')
ax2.axvline(x=1,linestyle='dashed')
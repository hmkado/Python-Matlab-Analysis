# -*- coding: utf-8 -*-
"""

@author: JR
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.physics.mechanics import *

#2.1
print('#2.1')
k=10*10**3
m=331
c=50
v=(65*1.60934*1000)/3600 #to m/s
Y=.05
zeta=c/(2*np.sqrt(k*m))
wb=v*2*np.pi
wn=np.sqrt(k/m)
r=wb/wn
X=Y*np.sqrt((1+(2*zeta*r)**2)/((1-r**2)**2+(2*zeta*r)**2))
theta1=np.arctan((2*zeta*wn*wb)/(wn**2-wb**2))
theta2=np.arctan((wn**2)/(wn**2-wb**2))
print('x(t) =',X,'cos(',wb,'*t','-',theta1,'-',theta2)

#2.2
print('\n#2.2')
n=100000
t=np.zeros(n+2)
tf=10
t[1]=0
delt=(tf-t[1])/n

x=np.zeros(n+2)
v=np.zeros(n+2)
x[1]=0
v[1]=0
for i in range(n):
    a=(c/m)*(Y*wb*np.cos(wb*t[i+1])-v[i+1])+(k/m)*(Y*np.sin(wb*t[i+1])-x[i+1])
    x[i+2]=x[i+1]+v[i+1]*delt
    v[i+2]=v[i+1]+a*delt
    t[i+2]=t[i+1]+delt
    
fig, (ax1, ax2) = plt.subplots(2)

ax1.plot(t, x)
ax1.set_title('Position [m]')
ax2.plot(t, v)
ax2.set_title('Velocity [m/s]')

#2.3
print('\n#2.3')
def xt(t):
    xp=X*np.cos(wb*t-theta1-theta2) #answer x(t)
    return xp
t=np.linspace(0,10,100)
plt.figure(2)
plt.title('displacement x (m) (linspace 100)')
plt.plot(t, xt(t))

t=np.linspace(0,10,100000)
plt.figure(3)
plt.title('displacement x (m) (linspace 100000)')
plt.plot(t, xt(t))
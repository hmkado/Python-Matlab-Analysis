# -*- coding: utf-8 -*-
"""

@author: JR
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.physics.mechanics import *

#3.1
print('#3.1')
dist=25
k=10*10**3
m=331
wn=np.sqrt(k/m)
v=(45*1.60934*1000)/3600 #to m/s
Y=.05
t1=dist/v #sec to reach to step
print(t1,'s to hit the step')
t=np.linspace(0,10,100)
d=v*t
def yy(t):
    y=np.piecewise(t, [t < 0, t >= 0], [0, Y])
    return y

f=k*yy(t-t1)
plt.figure(2)
plt.title('Force (N) input v. t')
plt.step(t,f)

#3.2
print('\n#3.2')

# settling ts=-ln(.03)/(zeta*wn) and ts=3, therefore c= ((-ln(.03)*2*sqrt(k*m)))/(3*wn)
ts=3
c=(-np.log(.03)*2*np.sqrt(k*m))/(ts*wn)
print('damping c=',c)

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
    a=(c/m)*(-v[i+1])+(k/m)*(yy(t[i+1]-t1)-x[i+1])
    x[i+2]=x[i+1]+v[i+1]*delt
    v[i+2]=v[i+1]+a*delt
    t[i+2]=t[i+1]+delt

fig, (ax1, ax2) = plt.subplots(2)

ax1.set_title('Position [m] v. time [s]')
ax1.plot(t, x)
ax1.plot(t, yy(t-t1))
ax1.axvline(3+t1, linestyle='dashed')
ax1.axhline(.05*(1.03), linestyle='dashed')
ax1.axhline(.05*(.93), linestyle='dashed')
ax2.plot(t, v)
ax2.set_title('Velocity [m/s]')

#3.3
print('\n#3.3')
t=np.linspace(0,10,10000)
zeta=c/(2*np.sqrt(k*m)) #zeta, c, wd from ts
if zeta < 1 : print('zeta = ',zeta,'less than 1, underdamped')
wd=wn*np.sqrt(1-zeta**2)
F0=k*Y
theta=np.arctan(zeta/np.sqrt(1-zeta**2))
x=((F0/k)-(F0/(k*np.sqrt(1-zeta**2)))*np.exp(-zeta*wn*(t-t1))*np.cos(wd*(t-t1)-theta))*np.heaviside(t-t1,1)
plt.figure(4)
plt.title('Analytical Position [m]')
plt.plot(t,x)

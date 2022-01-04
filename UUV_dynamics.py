"""
Created on Mon Jan  3 09:56:28 2022

@author: dgalli
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tikzplotlib import save as tikz_save

# Inputs
m = 40 # [kg]
Xud = -0.05 # [m]
r = 0.15 # [m]
Af = np.pi*r**2 # [m^2]
Cd = 0.05
rho = 1025.9 # [kg/m^3]

dt = 0.1

# Boundary conditions
t0 = 0   # [sec]
t1 = 600 # [sec]
q0 = 0   # [m]
q1 = 25  # [m]
v0 = 0   # [m/s]
v1 = 0   # [m/s]
a0 = 0   # [m-s2]
a1 = 0   # [m/s2]
j0 = 0
j1 = 0

DT = t1 - t0
h = q1 - q0

# Select the trajectory type
trajectory = "seventh"

# Init the lists
t = [t0]
q = [q0]
qd = [0]
qdd = [0]
Xuu = [0]
T = [0]

if trajectory == "linear": # q = a0 + a1(t - t0)
    method = "Linear Trajectory"
    a0 = q0
    a1 = h/DT
    for i in range(int(t1/dt)):
        q.append(a0 + a1*(t[i]-t0))
        u = h/DT
        ud = 0
        t.append(t[i] + dt)
        qd.append(u)
        qdd.append(ud)
        # Compute the drag
        Xuu.append(-0.5 * rho * Cd * Af * u * np.abs(u))
        # Compute the thrust
        T.append((m - Xud) * ud - Xuu[i])
        
    # Plot
    EqLab = "$q(t) = a_0 + a_1 \cdot (t-t_0)$"
    
if trajectory == "parabolic":
    method = "Parabolic Trajectory"
    tf = (t0 + t1)/2
    qf = (q0 + q1)/2
    
    a0 = q0
    a1 = v0
    a2 = (2/(DT**2))*(h - v0*DT)
    
    a3 = qf
    a4 = (2*h/DT) - v1
    a5 = (2/(DT**2))*(v1*DT - h)
    '''
    a0 = q0
    a1 = v0
    a2 = (4*h - DT*(3*v0+v1))/((2*DT)**2)
    a3 = (4*(q0+q1)+DT*(v0-v1))/8
    a4 = (4*h-DT*(v0+v1))/(2*DT)
    a5 = (-4*h+DT*(v0+3*v1))/((2*DT)**2)
    '''
    for i in range(int(t1/dt)):
        if t[i] < tf:
            q.append(a0 + a1*(t[i]-t0) + a2*((t[i]-t0)**2))
            u = a1 + 2*a2*(t[i]-t0)
            ud = 2*a2

        else:
            q.append(a3 + a4*(t[i]-tf) + a5*((t[i]-tf)**2))
            u = a4 + 2*a5*(t[i]-tf)
            ud = 2*a5
        t.append(t[i] + dt)
        qd.append(u)
        qdd.append(ud)
        # Compute the drag
        Xuu.append(-0.5 * rho * Cd * Af * u * np.abs(u))
        # Compute the thrust
        T.append((m - Xud) * ud - Xuu[i])
    EqLab = "$ q(t)=a_0 + a_1(t-t_0)+a_2(t-t_0)^2$  $t_0<t<t_f$ \n $ q(t)=a_3+a_4(t-t_f)+a_5(t-t_f)^2$  $t_f<t<t_1 $"
    
if trajectory == "cubic":
    method = "Cubic Trajectory"
    a0 = q0
    a1 = v0
    a2 = (3*h - DT*(2*v0+v1))/(DT**2)
    a3 = (-2*h+DT*(v0+v1))/(DT**3)
    
    for i in range(int(t1/dt)):
        q.append(a0 + a1*(t[i]-t0) + a2*((t[i]-t0)**2)+a3*((t[i]-t0)**3))
        u = a1 + 2*a2*(t[i]-t0) + 3*a3*((t[i]-t0)**2)
        ud = 2*a2 + 6*a3*(t[i]-t0)

        t.append(t[i] + dt)
        qd.append(u)
        qdd.append(ud)
        # Compute the drag
        Xuu.append(-0.5 * rho * Cd * Af * u * np.abs(u))
        # Compute the thrust
        T.append((m - Xud) * ud - Xuu[i])
        
    EqLab = "$ q(t)=a_0 + a_1(t-t_0)+a_2(t-t_0)^2+a_3(t-t_0)^3 $"
    
if trajectory == "fifth":
    method = "Fifth order Tajectory"
    a0 = q0
    a1 = v0
    a2 = 0.5*a0
    a3 = (1/(2*DT**3))*(20*h-(8*v1+12*v0)*DT-(3*a0-a1)*DT**2)
    a4 = (1/(2*DT**4))*(-30*h+(14*v1+16*v0)*DT-(3*a0-2*a1)*DT**2)
    a5 = (1/(2*DT**5))*(12*h-6*(v1+v0)*DT-(a1-a0)*DT**2)
    
    for i in range(int(t1/dt)):
        q.append(a0 + a1*(t[i]-t0) + a2*((t[i]-t0)**2)+a3*((t[i]-t0)**3)+a4*((t[i]-t0)**4)+a5*((t[i]-t0)**5))
        u = a1 + 2*a2*(t[i]-t0) + 3*a3*((t[i]-t0)**2) + 4*a4*((t[i]-t0)**3) + 5*a5*((t[i]-t0)**4)
        ud = 2*a2 + 6*a3*(t[i]-t0) + 12*a4*((t[i]-t0)**2) + 20*a5*((t[i]-t0)**3)

        t.append(t[i] + dt)
        qd.append(u)
        qdd.append(ud)
        # Compute the drag
        Xuu.append(-0.5 * rho * Cd * Af * u * np.abs(u))
        # Compute the thrust
        T.append((m - Xud) * ud - Xuu[i])
        
    EqLab = "$ q(t)=a_0 + a_1(t-t_0)+a_2(t-t_0)^2+a_3(t-t_0)^3+a_4(t-t_0)^4+a_5(t-t_0)^5$"

''' It does not work properly
if trajectory == "seventh":
    method = "Seventh order Trajectory"
    a0 = q0
    a1 = v0
    a2 = 0.5*a0
    a3 = j0/6
    a4 = ( 210*h - DT*((30*a0 - 15*a1)*DT + (4*j0 +   j1)*DT**2 + 120*v0 +  90*v1))/(6*DT**4)
    a5 = (-168*h + DT*((20*a0 - 14*a1)*DT + (2*j0 +   j1)*DT**2 +  90*v0 +  78*v1))/(2*DT**5)
    a6 = (  42*h - DT*((45*a0 - 39*a1)*DT + (4*j0 + 3*j1)*DT**2 + 216*v0 + 204*v1))/(6*DT**6)
    a7 = ( 210*h - DT*((12*a0 - 12*a1)*DT + (  j0 +   j1)*DT**2 +  60*v0 +  60*v1))/(6*DT**7)
    
    for i in range(int(t1/dt)):
        q.append(a0 + a1*(t[i]-t0) + a2*((t[i]-t0)**2) + a3*((t[i]-t0)**3) + a4*((t[i]-t0)**4) + a5*((t[i]-t0)**5) + a6*((t[i]-t0)**6) + a7*((t[i]-t0)**7))
        u = a1 + 2*a2*(t[i]-t0) + 3*a3*((t[i]-t0)**2) + 4*a4*((t[i]-t0)**3) + 5*a5*((t[i]-t0)**4) + 6*a6*((t[i]-t0)**5) + 7*a7*((t[i]-t0)**6)
        ud = 2*a2 + 6*a3*(t[i]-t0) + 12*a4*((t[i]-t0)**2) + 20*a5*((t[i]-t0)**3) + 30*a6*((t[i]-t0)**4) + 42*a7*((t[i]-t0)**5)

        t.append(t[i] + dt)
        qd.append(u)
        qdd.append(ud)
        # Compute the drag
        Xuu.append(-0.5 * rho * Cd * Af * u * np.abs(u))
        # Compute the thrust
        T.append((m - Xud) * ud - Xuu[i])
        
    EqLab = "$ q(t)=a_0 + a_1(t-t_0)+a_2(t-t_0)^2+a_3(t-t_0)^3+a_4(t-t_0)^4+a_5(t-t_0)^5+a_6(t-t_0)^6+a_7(t-t_0)^7$"
'''

# Plot the outcomes

fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize = (10, 15))
fig.suptitle(method)
ax1.plot(t, q, color = "red", label = EqLab)
ax2.plot(t, qd, color = "blue", label = "$\dot{q}$ $[m/s]$")
ax3.plot(t, qdd, color = "green", label = "$\ddot{q}$ $[m/s^2]$")
ax4.plot(t, T, color = "black",  label = "thrust $T$ $[N]$")
ax4.plot(t, Xuu, color = "orange", label = "drag $X_{u|u|}$ $[N]$")
plt.xlabel("time $t [s]$")

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

ax1.set_title('Position', fontstyle='italic')
ax2.set_title('Velocity', fontstyle='italic')
ax3.set_title('Acceleration', fontstyle='italic')
ax4.set_title('Thrust and Drag forces', fontstyle='italic')


ax1.set_xlim(t0, t1)
ax2.set_xlim(t0, t1)
ax3.set_xlim(t0, t1)
ax4.set_xlim(t0, t1)
ax1.set_ylim(q0,q1)
plt.show()
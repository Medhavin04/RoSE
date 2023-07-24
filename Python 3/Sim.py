import numpy as np
from numpy import pi
from numpy.linalg import inv
from numpy import matmul
from numpy import identity as I

from scipy import signal
import matplotlib.pyplot as plt

import csv

P0 = 101325

x = y = vx = vy = ax = ay = t = 0
dt = 0.05
g = 9.8

launch_angle = 15*pi/180

burn_increase = 10
burn_steady = 10
burn_decrease = 10

max_acceleration = 20

X = np.array([])
Y = Z = np.array([])

def alt(P):
    
    return round(4947.18853894153 * (pow(P0,0.190255) - pow(P,0.190255)),2)

def motor(t):
    
    if t > burn_increase + burn_steady + burn_decrease:
        return 0
    if t < burn_increase:
        return max_acceleration*(t)/burn_increase
    if burn_increase <= t <= burn_increase + burn_steady:
        return max_acceleration
    if burn_increase + burn_steady < t <= burn_increase + burn_steady + burn_decrease:
        return max_acceleration*( + burn_increase + burn_steady + burn_decrease - t)/burn_decrease

v1 = v2 = 0

with open('data.csv','w') as file:
    csvwriter = csv.writer(file)
    
    while y >= 0:

        if t <  + burn_increase + burn_steady + burn_decrease:
            ax = motor(t)*np.sin(launch_angle)
            ay = motor(t)*np.cos(launch_angle)
        else:
            ax = 0
            ay = -g
            
        x += (vx+ax/2*dt)*dt
        y += (vy+ay/2*dt)*dt
        
        vx += ax*dt
        vy += ay*dt
        
        if ax <= 0 and t > burn_increase + burn_steady + burn_decrease:
            pitch = pi/2 - np.arctan(vy/vx)
        else:
            pitch = launch_angle
        
        P = pow(pow(P0,0.190255) - y/4947.18853894153,1/0.190255)
        a1 = ax*np.cos(pitch) - ay*np.sin(pitch)
        a2 = ax*np.sin(pitch) + ay*np.cos(pitch)
        a3 = 0
        
        v1 += a1*dt
        v2 += a2*dt
        
        data = [P+np.random.normal(0,1),a1 + np.random.normal(0,0.5),a2 + np.random.normal(0,0.5),a3 + np.random.normal(0,0.5),P,a1,a2,a3,pitch,max_acceleration,burn_increase,burn_steady,burn_decrease]
        
        csvwriter.writerow(data)
        
#         X = np.append(X,t)
#         Y = np.append(Y,pitch*180/pi)
        
        t += dt

file.close()

# for i in range(len(X)):
#         plt.plot(X[i],Y[i], marker="o", markersize=1, markeredgecolor="red", markerfacecolor="blue")
# 
# plt.show()
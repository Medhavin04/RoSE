# CHANGES FROM v1:
# 1. REMOVED LPF AFTER KF
# 2. INTEGRATED ALL VALUES (BOTH MEASURED & DERIVED INTO ONE MATRIX)
# 3. ADDED PITCH ESTIMATION (DOESN'T WORK WELL WHEN ROCKET IS NEAR VERITCAL (10 DEGRESS) OR SPEED IS LOW)
# 4. ADDED GRAPHS
# 5. ADDED SIMULATOR TO GENERATE SAMPLE NOISY DATA (Sim.py)

import numpy as np
from numpy.linalg import inv
from numpy import matmul
from numpy import identity as I

from scipy import signal
import matplotlib.pyplot as plt
from numpy import pi

# KF ----------------------------------------------------------------------------------------------------------------------------------------------

def alt(P):
    
    return round(4947.18853894153 * (pow(Data[0][0],0.190255) - pow(P,0.190255)),2)

def motor(t):         # THRUST CURVE
    
    if t < Data[0][10]:
        return Data[0][9]*t
    if Data[0][10] <= t <= Data[0][10]+Data[0][11]:
        return Data[0][9]
    if Data[0][10]+Data[0][11] < t < Data[0][10]+Data[0][11]+Data[0][12]:
        return Data[0][9]*(Data[0][10]+Data[0][11]+Data[0][12] - t)/Data[0][12]
    else:
        return 0

counter = 0 

# COPY DATA FROM CSV TO Data[]

Data = np.loadtxt('data.csv', delimiter=',')

# KALMAN FILTER MATRICES

Xm = np.full((9,1),0.)    # MEASUREMENTS MATRIX
X = np.full((9,1),0.)     # STATE MATRIX
W = np.full((9,1),0.)     # STATE PREDITCION NOISE MATRIX (0 SEEMS TO WORK FINE)
Q = np.full((9,9),0.)     # STATE UNCERTAINITY PREDICTION NOISE MATRIX
R = np.full((9,9),0.)     # SENSOR NOISE COVARIANCE MATRIX
P = np.full((9,9),0.)     # STATE UNCERTAINITY MATRIX
Xp = X                    # PREDICTED STATE MATRIX
Pp = P                    # PREDICTED STATE UNCERTAINITY MATRIX
K = np.full((9,9),0.)     # KALMAN GAIN MATRIX
U = np.full((9,1),0.)     # CONTROL MATRIX

# INITIALIZING STATE MATRIX WITH FIRST VALUES OF DATA (MAKE SURE INITIAL STATE IS STATIONARY)

X[0][0] = Data[0][0]
X[1][0] = Data[0][1]
X[2][0] = Data[0][2]
X[3][0] = Data[0][3]
X[4][0] = X[5][0] = X[6][0] = X[7][0] = X[8][0] = 0

# SENSOR NOISE MATRIX VALUES

R[0][0] = 2
R[1][1] = 0.01
R[2][2] = 0.01
R[3][3] = 0.01
R[4][4] = 0.01
R[5][5] = 5
R[6][6] = 800
R[7][7] = 0.2
R[8][8] = 1

# SETTING STATE PREDICTION UNCERTAINITY

Q[0][0] = 0.4
Q[1][1] = 0.005
Q[2][2] = 0.005
Q[3][3] = 0.005
Q[4][4] = 0.1
Q[5][5] = 2
Q[6][6] = 50
Q[7][7] = 0.01
Q[8][8] = 2

P = Q

dt = 0.05         # FIRST TIME STEP
g = 9.8           # ACCELERATION DUE TO GRAVITY

A = Ar = As = np.array([])

c = 1

# FILTER

for counter in range(len(Data)-1):
    
    # UPDATE CONTROL MATRIX
    U[5][0] = X[6][0]*dt
    if counter*dt > 1:
        U[7][0] = (X[5][0]*X[2][0]*dt - X[8][0]*U[5][0])/(X[8][0]*X[8][0]*np.sin(X[7][0]))
        
    U[0][0] = -(X[5][0]+X[6][0]/2*dt)*dt/(4947.18853894153*0.190255*pow(X[0][0],0.190255-1))
    U[1][0] = g*np.cos(X[7][0])*U[7][0]
    U[2][0] = motor((counter+1)*dt) - motor(counter*dt) + g*np.sin(X[7][0])*U[7][0]
    U[3][0] = 0
    U[4][0] = (X[5][0]+X[6][0]/2*dt)*dt
    U[6][0] = U[2][0]*np.cos(X[7][0])

    # UPDATE STATE PREDICTION & STATE UNCERTAINITY PREDICTION

    Xp = X + U
    
    Pp = P + Q
    
    # UPDATE KALMAN GAIN MATRIX
    
    K = matmul(Pp,inv(Pp+R))
    
    # UPDATE MEASUREMENTS MATRIX
    
    Xm[0][0] = Data[counter][0]
    Xm[1][0] = Data[counter][1]
    Xm[2][0] = Data[counter][2]
    Xm[3][0] = Data[counter][3]
    
    last_alt = Xm[4][0]
    Xm[4][0] = alt(Xm[0][0])
    last_vel = Xm[5][0]
    Xm[5][0] = c*(Xm[4][0] - last_alt)/dt + (1-c)*Xm[8][0]*np.cos(X[7][0])
    
    if abs(np.cos(Xm[7][0])) > 0.01:
        Xm[6][0] = Xm[2][0]/np.cos(Xm[7][0])
    
    Xm[8][0] += Xm[2][0]*dt
    
    if abs(Xm[5][0]) <= abs(Xm[8][0]):
        Xm[7][0] = np.arccos(Xm[5][0]/Xm[8][0])
    
    # UPDATE STATE & STATE UNCERTAINITY
    
    X = Xp + matmul(K,(Xm - Xp))
    
    P = matmul((I(len(K)) - K),Pp)
    
#     Ar = np.append(Ar,Xm[7][0]*180/pi)
#     A = np.append(A,X[7][0]*180/pi)
#     As = np.append(As,Data[counter][8]*180/pi)
    
    # MOVE TO NEXT READINGS

    counter += 1
    
    # UPDATE dt
    
    dt = 0.05
    
    # REPEAT TILL EOF

# GRAPHS --------------------------------------------------------------------------------------

# tlims = [0,0.05*len(Data)]        # in seconds
# t = np.linspace(tlims[0],tlims[1],len(Data)-1)
#     
# # Plot the signal
# plt.figure()
# plt.plot(t,A)
# plt.plot(t,Ar)
# plt.plot(t,As)
# plt.ylabel("$y(t)$")
# plt.xlim([0,dt*len(Data)]);
# 
# # # Generate Fourier transform
# # yfilthat = np.fft.fft(yfilt)
# # fcycles = np.fft.fftfreq(len(t),d=1.0/samplingFreq)
# # 
# # plt.figure()
# # plt.plot(fcycles,np.absolute(yhat));
# # plt.plot(fcycles,np.absolute(yfilthat));
# # plt.xlim([-7.5,7.5]);
# # plt.xlabel("$\omega$ (cycles/s)");
# # plt.ylabel("$|\hat{y}|$");
# 
# plt.show()
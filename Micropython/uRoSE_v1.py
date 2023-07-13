from ulab import numpy as np

def I(n):
    
    identity = np.full((n,n),0.)
    
    for i in range(n):
        
        identity[i][i] = 1.
        
    return identity

# KF + (OPTIONAL) LPF ----------------------------------------------------------------------------------------------------------------------------------------------

def altitude(P,T):    # RETURNS (ALT FROM SEA LEVEL - ALT OF GROUND FROM SEA LEVEL)
    
    return 1/0.0065 * (T + 273.15) * (pow(101325/P,1/5.257) - 1) - 1/0.0065 * (Data[0][8] + 273.15) * (pow(101325/Data[0][9],1/5.257) - 1)

def motor(t):         # RATE OF CHANGE OF ACCELERATION FROM THRUST CURVE (FOR CONTROL MATRIX)
    
    if 60750 <= t <= 61369:
        
        return 75/0.619
    
    if 61369 < t < 61992:
        
        return 0
    
    if 61992 <= t <= 62665:
        
        return -75/0.673
    
    else:
        
        return 0

counter = 0 

# COPY DATA FROM CSV TO Data[]

Data = np.loadtxt('launch_3_data.csv', delimiter=',')

# KALMAN FILTER MATRICES

Xm = np.full((5,1),0.)    # MEASUREMENTS MATRIX
X = np.full((5,1),0.)     # STATE MATRIX
W = np.full((5,1),0.)     # STATE PREDITCION NOISE MATRIX (0 SEEMS TO WORK FINE)
Q = np.full((5,5),0.)     # STATE UNCERTAINITY PREDICTION NOISE MATRIX
R = np.full((5,5),0.)     # SENSOR NOISE COVARIANCE MATRIX
P = np.full((5,5),0.)     # STATE UNCERTAINITY MATRIX
Xp = X                    # PREDICTED STATE MATRIX
Pp = P                    # PREDICTED STATE UNCERTAINITY MATRIX
K = np.full((5,5),0.)     # KALMAN GAIN MATRIX
U = np.full((5,1),0.)     # CONTROL MATRIX

# LPF AND BASIC SENSOR FUSION MATRICES

Xd1 = np.full((3,1),0.)   # DERIVED VALUES MATRIX (FROM BAROMETER)
Xd2 = np.full((3,1),0.)   # DERIVED VALUES MATRIX (FROM ACCELEROMETER)

# DERIVED VALUES - ALTITUDE, VELOCITY, ACCELERATION

Xf = np.full((3,1),0.)    # FUSED VALUES MATRIX
Xl = np.full((3,1),0.)    # RESULT OF LPF ON Xf
G1 = np.full((3,3),0.)    # WEIGHTAGE OF BAROMETER vs ACCELEROMETER
G2 = np.full((3,3),0.)    # LPF WEIGHTAGE

# INITIALIZING STATE MATRIX WITH FIRST VALUES OF DATA (MAKE SURE INITIAL STATE IS STATIONARY)

X[0][0] = Data[0][9]
X[1][0] = Data[0][8]
X[2][0] = Data[0][2]
X[3][0] = Data[0][3]
X[4][0] = Data[0][4]

# SENSOR NOISE MATRIX VALUES

R[0][0] = 4
R[1][1] = 0.09
R[2][2] = 0.01
R[3][3] = 0.01
R[4][4] = 0.01

# INITIALZING STATE UNCERTAINITY

P = R

# SETTING STATE PREDICTION UNCERTAINITY

Q[0][0] = 0.4
Q[1][1] = 0.09
Q[2][2] = 0.015
Q[3][3] = 0.015
Q[4][4] = 0.015

dt = 0.05         # FIRST TIME STEP
g = 9.8           # ACCELERATION DUE TO GRAVITY

# SETTING BARO vs ACCELEROMETER WEIGHTAGES (HIGHER THE VALUE, MORE THE TRUST IN BAROMETER)

G1[0][0] = 1
G1[1][1] = 0.75
G1[2][2] = 0

# SETTINGS LPF WIGHTAGES (MORE THE VALUE, FASTER & NOISIER THE RESPONSE)
# IF VALUE IS 1, THERE IS EFFECTIVELY NO LPF. I DON'T THINK WE NEED ONE FOR RIGHT NOW, BUT WE MIGHT FOR PITCH ESTIMATION BECAUSE NOISE WILL REALLY MESS THAT UP

G2[0][0] = 1
G2[1][1] = 1
G2[2][2] = 1

# FILTER

for counter in range(len(Data)-1):
    
    # UPDATE CONTROL MATRIX
    
    U[0][0] = - (0.0341705*pow(X[0][0],1+1/5.257)*(2 - pow(101325/X[0][0],1/5.257))*(Xl[1][0] + Xl[2][0]/2*dt)*dt)/((X[1][0] + 273.15)*pow(101325,1/5.257))
    U[1][0] = - 0.0001*(Xl[1][0] + Xl[2][0]*dt)*dt
    U[2][0] = motor(Data[0][0])*dt

    # UPDATE STATE PREDICTION & STATE UNCERTAINITY PREDICTION

    Xp = X + U
    
    Pp = P + Q
    
    # UPDATE KALMAN GAIN MATRIX
    
    K = np.dot(Pp,np.linalg.inv(Pp + R))
    
    # UPDATE MEASUREMENTS MATRIX
    
    Xm[0][0] = Data[counter][9]
    Xm[1][0] = Data[counter][8]
    Xm[2][0] = -Data[counter][2] - g
    Xm[3][0] = Data[counter][3]
    Xm[4][0] = Data[counter][4]
    
    # UPDATE STATE & STATE UNCERTAINITY
    
    X = Xp + np.dot(K,(Xm - Xp))
    
    P = np.dot((I(len(K)) - K),Pp)
    
    # UPDATE DERIVED VALUES FROM BAROMETER
    
    Xd1[0][0] = altitude(X[0][0],X[1][0])     # ALTITUDE CALCULATED USING KF OUTPUT OF PRESSURE AND TEMP
    Xd1[1][0] = (Xd1[0][0] - Xl[0][0])/dt     # VELOCTY = d(ALTITUDE)/dt
    Xd1[2][0] = (Xd1[1][0] - Xl[1][0])/dt      # ACCELERATION = d(VELOCTY)/dt
    
    # UPDATE DERIVED VALUES FROM ACCELEROMETER (DEAD RECKONING)
    
    Xd2[2][0] = X[2][0]                       # ACCELERATION = OUTPUT FROM KF
    Xd2[1][0] += Xd2[2][0]*dt                 # VELOCITY = PREVIOUS VELOCITY + ACCELERATION*dt
    Xd2[0][0] += Xd2[1][0]*dt                 # ALTITUDE = PREVIOUS ALTITUDE + VELOCITY*dt
    
    # FUSE Xd1 & Xd2
    
    Xf = np.dot(G1,Xd1) + np.dot((I(len(G1)) - G1),Xd2)
    
    # APPLY LPF ON Xf

    Xl = np.dot(G2,Xf) + np.dot((I(len(G2)) - G2),Xl)
    
    # FILTERED vs RAW ALTITUDE - print(Xl[0][0],altitude(Data[counter][9],Data[counter][8]))

    # FILTERED vs RAW VELOCITY - print(Xl[1][0],(altitude(Data[counter][9],Data[counter][8]) - altitude(Data[counter-1][9],Data[counter-1][8]))/dt)

    # FILTERED vs RAW ACCELERATION - print(Xl[2][0],-float(Data[counter][2])-g)
    
    
    
    # MOVE TO NEXT READINGS

    counter += 1
    
    # UPDATE dt
    
    dt = (Data[counter][0] - Data[counter - 1][0])/1000
    
    # REPEAT TILL EOF


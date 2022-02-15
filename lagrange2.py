from scipy.optimize import root_scalar
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

G = 6.67408e-11
m1 = 1.9885e30 # kg
m2 = 5.972e24 # kg
R = 1.49598e8 # km
w = 1.992e-7

def LN(r):
    f1 = -G*m1/((r)*abs(r))
    f2 = -G*m2/((r-R)*abs(r-R))
    ang_v = G*r*(m1+m2)/(R**3)
    return ang_v + f1 + f2


L2 = root_scalar(LN, bracket=[1.50e8, 1.52e8]).root

def r1(x,y,z):
    return np.sqrt(x**2 + y**2 + z**2)


def r2(x,y,z):
    return np.sqrt((x-R)**2 + y**2 + z**2)


def EOM(t, X):
    Xdot = np.zeros(len(X))
    Xdot[:3] = X[3:6]
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    R1 = r1(x,y,z)
    R2 = r2(x,y,z)
    xddot = x*w**2 + 2*w*ydot - G*m1*x/(R1**3) - G*m2*(x-R)/(R2**3)
    yddot = y*w**2 - 2*w*xdot - G*m1*y/(R1**3) - G*m2*y/(R2**3)
    zddot = - G*m1*z/(R1**3) - G*m2*z/(R2**3)
    Xdot[3] = xddot
    Xdot[4] = yddot
    Xdot[5] = zddot
    return Xdot


def main():
    CONDITION_0 = np.array([L2 + 100/R, 0, 0, 0, 0, 0])
    EOMs = EOM
    

if __name__ == '__main__':
    main()
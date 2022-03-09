from tabnanny import check
from scipy.optimize import root_scalar
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

G = 6.67408e-11
m1 = 5.972e24 # 1.9885e30 # kg
m2 = 7.35e22 #5.972e24 # kg
R = 3.844e8 #1.49598e11 # m
radial_velocity = np.sqrt(G*(m1+m2)/R**3)

def LN(r):
    f1 = -G*m1/((r)*abs(r))
    f2 = -G*m2/((r-R)*abs(r-R))
    ang_v = G*r*(m1+m2)/(R**3)
    return ang_v + f1 + f2


#L2 = root_scalar(LN, bracket=[1.50e11, 1.52e11]).root
L2 = root_scalar(LN, bracket=[0.448e9, 0.45e9]).root

def solve_lde(f, variable):
    ...


def r1(x,y,z):
    return np.sqrt(x**2 + y**2 + z**2)


def r2(x,y,z):
    return np.sqrt((x-R)**2 + y**2 + z**2)


def nonlinear_EOM(t, X):
    Xdot = np.zeros(len(X))
    Xdot[:3] = X[3:6]
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    R1 = r1(x,y,z)
    R2 = r2(x,y,z)
    xddot = x*radial_velocity**2 + 2*radial_velocity*ydot - G*m1*x/(R1**3) \
        - G*m2*(x-R)/(R2**3)
    yddot = y*radial_velocity**2 - 2*radial_velocity*xdot - G*m1*y/(R1**3) \
        - G*m2*y/(R2**3)
    zddot = - G*m1*z/(R1**3) - G*m2*z/(R2**3)
    Xdot[3] = xddot
    Xdot[4] = yddot
    Xdot[5] = zddot
    return Xdot


def init_vy(lpoint, distance):
    x= lpoint + distance
    a = x*G*(m1+m2)/(R**3) - x*G*m1/(x**3) - G*m2*(R-x)/(x**3)
    print(a)
    return np.sqrt(a*distance)

def main2():
    fig, ax = plt.subplots()
    x = np.arange(-2*R, 2*R, 10000)
    y = LN(x)
    ax.plot(x,y)
    ax.set_ylim(-0.01, 0.01)
    plt.show()

def main():
    ixdot = 0
    iydot = 0.0014777#init_vy(L2, 1000) # ms-1
    print(init_vy(L2, 100))
    CONDITION_0 = np.array([L2 - 100, 0, 0, ixdot, iydot, 0])
    EOMs = nonlinear_EOM
    traj = solve_ivp(EOMs, [0,28*24*3600], CONDITION_0, atol=0.000001, rtol=3e-14)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(traj.y[0,:], traj.y[1,:], traj.y[2,:], 'b')
    ax.plot(L2, 0, 0, 'k+')
    ax.plot(L2-100, 0, 0, 'r+')

    bound = 500
    ax.axes.set_xlim3d(left=L2 - bound, right=L2 + bound)
    ax.axes.set_ylim3d(bottom=-bound, top=bound)
    ax.axes.set_zlim3d(bottom=-bound, top=bound)

    plt.show()

if __name__ == '__main__':
    main()
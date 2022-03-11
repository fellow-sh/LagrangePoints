from scipy.optimize import root_scalar
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

G = 6.67408e-11
m1 = 1.9885e30 # kg
m2 = 5.972e24 # kg
R = 1.49598e11 # m
radial_velocity = np.sqrt(G*(m1+m2)/R**3)

def LN(r):
    f1 = -G*m1/((r)*abs(r))
    f2 = -G*m2/((r-R)*abs(r-R))
    ang_v = G*r*(m1+m2)/(R**3)
    return ang_v + f1 + f2


L2 = root_scalar(LN, bracket=[1.50e11, 1.52e11]).root
#L2 = root_scalar(LN, bracket=[0.448e9, 0.45e9]).root


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


def init_circular_v(lpoint, position: float):
    x= lpoint + position
    R1 = r1(x, 0, 0)
    R2 = r2(x, 0, 0)
    xddot = x*radial_velocity**2 - G*m1*x/(R1**3) \
        - G*m2*(x-R)/(R2**3)
    print(f'a_r: {xddot}')
    ydot = 2.2*np.sqrt(xddot*position)
    return ydot


def main2():
    fig, ax = plt.subplots()
    x = np.arange(-2*R, 2*R, 10000)
    y = LN(x)
    ax.plot(x,y)
    ax.set_ylim(-0.01, 0.01)
    plt.show()


def main():
    time_u = 24*3600 # time unit (days)
    duration = ... # days
    distance = 1000 # m

    def plot1() -> None:
        duration = 50 # days
        YDOT_0 = 0 # ms-1
        CONDITION_0 = np.array([L2 - distance, 0, 0, 0, YDOT_0, 0])
        traj = solve_ivp(nonlinear_EOM, [0,time_u*duration], CONDITION_0, atol=1e-6, rtol=3e-14)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot(traj.y[0,:], traj.y[1,:], traj.y[2,:])
        ax.plot(L2, 0, 0, 'k+')
        ax.plot(L2-distance, 0, 0, 'r+')
        ax.text(L2-distance, 0, 0, 'P₀')

        #ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        bound = 3.5*distance
        ax.axes.set_xlim3d(left=L2 - bound, right=L2 + bound)
        ax.axes.set_ylim3d(bottom=-bound, top=bound)
        ax.axes.set_zlim3d(bottom=-bound, top=bound)
        ax.view_init(20, -140)

        plt.show()


    def plot2() -> None:
        ...


    plot1()

    iydot = init_circular_v(L2, distance)
    #iydot = .00000130549*distance
    iydot = 0
    #print(f'test_iydot: {iydot}', f'c / test: {iydot/init_circular_v(L2, distance)}')
    # precision in the initial conditions is crucial
    #iydot = 0.147763 # ms-1

    plt.show()


if __name__ == '__main__':
    main()

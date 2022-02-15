from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
import numpy as np

G = 6.67408e-11
m1 = 1.9885e30 # kg
m2 = 5.972e24 # kg
R = 1.49598e8 # km

def LN(r):
    #TODO: derive this equation mathematically (checked!)
    f1 = -G*m1/((r)*abs(r))
    f2 = -G*m2/((r-R)*abs(r-R))
    ang_v = G*r*(m1+m2)/(R**3)
    #ang_v = G*((m1+m2)*(x3+x1)-(m2*l))/(l**3)
    return ang_v + f1 + f2

l1 = root_scalar(LN, bracket=[1.48e8, 1.49e8]).root
l2 = root_scalar(LN, bracket=[1.50e8, 1.52e8]).root
l3 = root_scalar(LN, bracket=[-1.52e8, -1.49e8]).root

print('L1:', l1)
print('L2:', l2)
print('L3:', l3)

domain = np.arange(-2*R, 2*R, 100)
range_ = LN(domain)

fig, ax = plt.subplots(1, 2, figsize=(8.5,4), constrained_layout=True)
ax[0].plot(domain, range_)
ax[0].grid()
ax[0].set_ylim(-0.25e5, 0.25e5)
ax[0].set_xlim(-abs(l3)*1.2, abs(l3*1.2))
#ax[0].set_title('Net radial acceleration in the Sun-Earth system')
ax[0].plot(0, 0, 'oy')
ax[0].plot(R, 0, 'og')
ax[0].plot(l3, 0, '.k')
ax[0].annotate('Sun',
               xy=(0,0), xycoords='data',
               xytext=(0,15), textcoords='offset points',
               horizontalalignment='center', verticalalignment='top')
ax[0].annotate('Earth',
               xy=(R,0), xycoords='data',
               xytext=(-2,15), textcoords='offset points',
               horizontalalignment='right', verticalalignment='top')
ax[0].annotate('L3',
               xy=(l3,0), xycoords='data',
               xytext=(0,15), textcoords='offset points',
               horizontalalignment='center', verticalalignment='top')
#ax[0].set_xlabel(u'distance to m\u2081 (km)')
#ax[0].set_ylabel(u'radial acceleration (m s\u207B\u00B2)')

ax[1].plot(domain, range_)
ax[1].grid()
ax[1].set_xlim(l1/1.01, l2*1.01)
ax[1].set_ylim(-5000, 5000)
#ax[1].set_title('Net radial acceleration around Earth')
ax[1].plot(l1, 0, '.k')
ax[1].plot(l2, 0, '.k')
ax[1].plot(R, 0, 'og')
ax[1].annotate('L1',
               xy=(l1,0), xycoords='data',
               xytext=(0,15), textcoords='offset points',
               horizontalalignment='center', verticalalignment='top')
ax[1].annotate('L2',
               xy=(l2,0), xycoords='data',
               xytext=(0,15), textcoords='offset points',
               horizontalalignment='center', verticalalignment='top')
ax[1].annotate('Earth',
               xy=(R,0), xycoords='data',
               xytext=(0,15), textcoords='offset points',
               horizontalalignment='center', verticalalignment='top')
#ax[1].set_xlabel(u'distance to m\u2081 (km)')
#ax[1].set_ylabel(u'radial acceleration (m s\u207B\u00B2)')
plt.show()

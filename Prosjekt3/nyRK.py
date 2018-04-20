import numpy as np
import matplotlib.pyplot as plt

Er = 6371e3         #Earth radius [m]
Leo = 322e3         #Low earth orbit [m]
Geo = 35680e3       #Geostationary earth orbit [m]
r_leo = Leo+Er         #Distance from center of earth to satellite [m]
r_geo = Geo + Er
p = 120*60          #Orbital period [s]
G_const = 6.673e-11 #Gravitational constant [N (m/kg)^2]
M = 5.9726e24       #Earth mass [kg]
m = 720             #Satellite mass [kg]
v0_LEO = np.sqrt(G_const * M / r_leo)  # m/s
v0_GEO = np.sqrt(G_const * M / r_leo)
w0 = np.array([r_leo, 0, 0, v0_LEO])

t_max = 10000
dt = 1


def Runge_Kutta(f, t_max, dt, w0):

    timelist = np.arange(0, t_max + dt, dt)
    positionlist = np.zeros((len(timelist), 4))

    w = w0
    positionlist[0] = w

    for i in range(1, (len(timelist))):
        s1 = f(w)
        s2 = f(w + dt / 2 * s1)
        s3 = f(w + dt / 2 * s2)
        s4 = f(w + dt * s3)

        w = w + dt / 6 * (s1 + 2 * s2 + 2 * s3 + s4)
        positionlist[i] = w

    # Return x-values, y-values, vx-values, vy-values as separate vectors
    return positionlist[:, 0], positionlist[:, 2], positionlist[:, 1], positionlist[:, 3]

 #Derivative function (equation set) with SI units for trajectory
def f_SI(w):
    r = np.sqrt(w[0] ** 2 + w[2] ** 2)
    return np.array([w[1],
                     - G_const * M * w[0] / r ** 3,
                     w[3],
                     - G_const * M * w[2] / r ** 3])


# Runge-Kutta for single timestep and starting state
x, y, vx, vy = Runge_Kutta(f_SI, t_max, dt, w0)


plt.figure()
plt.title('LEO orbit')
plt.plot(x, y, 'ro', label = r'LEO satellite ($\epsilon=0.0$)',color = 'purple', linewidth = 0.1)
#plt.plot(GSx, GSy, 'g', label = r'Geo satellite ($\epsilon=0.0$)',color = 'orange', linewidth = 0.5)
plt.plot(0,0,'ro', label='Earth', color = 'blue')
plt.legend()
plt.xlabel(r"$x$ [r]")
plt.ylabel(r"$y$ [r]")
plt.axes().set_aspect('equal','datalim')
plt.grid()
plt.show()
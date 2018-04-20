import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Set common figure parameters:
newparams = {
    'figure.figsize': (16, 5), 'axes.grid': True,
    'lines.linewidth': 1.5, 'font.size': 19, 'lines.markersize': 10,
    'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral'}
plt.rcParams.update(newparams)

Er = 6371e3  # Earth radius [m]
Leo = 322e3  # Low earth orbit [m]
Geo = 35680e3  # Geostationary earth orbit [m]
r_leo = Leo + Er  # Distance from center of earth to satellite [m]
r_geo = Geo + Er

#p = 120 * 60  # Orbital period [s]
G_const = 6.673e-11  # Gravitational constant [N (m/kg)^2]
M = 5.9726e24  # Earth mass [kg]
m = 720  # Satellite mass [kg]
v0_LEO = np.sqrt(G_const * M / r_leo)  # m/s
v0_GEO = np.sqrt(G_const * M / r_leo)
w0 = np.array([r_leo, 0, 0, v0_LEO])
# C_leo = (G_const * M * p **2 /r_leo**3)
# C_geo = (G_const * M * p **2 /r_geo**3)

# print(C)

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

    # Derivative function (equation set) with SI units for trajectory


def f(w):
    r = np.sqrt(w[0] ** 2 + w[2] ** 2)
    return np.array([w[1],
                     - G_const * M * w[0] / r ** 3,
                     w[3],
                     - G_const * M * w[2] / r ** 3])




'''''''''
    plt.figure()
    plt.title('Energy as a function of time')
    plt.plot(time_vector, E_4RK)
    plt.ylabel("Energy [J]")
    plt.xlabel("Time [yr]")
    plt.grid()
    plt.show()

    plt.figure()
    plt.title('Potenial energy as a function of time')
    plt.plot(time_vector, P_4RK)
    plt.ylabel("Energy [J]")
    plt.xlabel("Time [yr]")
    plt.grid()
    plt.show()

    plt.figure()
    plt.title('Kinetic energy as a function of time')
    plt.plot(time_vector, K_4RK)
    plt.ylabel("Energy [J]")
    plt.xlabel("Time [yr]")
    plt.grid()
    plt.show()

    plt.figure()
    plt.title('Speed as a function of time')
    plt.plot(time_vector, total_velocity)
    plt.ylabel("Speed [AU/yr]")
    plt.xlabel("Time [yr]")

    plt.grid()
    plt.show()
    '''''''''''


def plotLEOtrajectory():
    # LEO sattelite trajectory for fixed earth
    x, y, vx, vy = Runge_Kutta(f, t_max, dt, w0)

    plt.figure()
    plt.title('LEO and GEO orbits')
    plt.plot(x, y, 'ro', label=r'LEO satellite ($\epsilon=0.0$)')
    # plt.plot(GSx, GSy, 'g', label = r'Geo satellite ($\epsilon=0.0$)',color = 'orange', linewidth = 0.5)
    plt.plot(0, 0, 'ro', label='Earth', color='blue')
    plt.legend(loc=1)
    plt.xlabel(r"$x$ [r]", size = 24)
    plt.ylabel(r"$y$ [r]", size = 24)
    plt.axes().set_aspect('equal', 'datalim')
    plt.grid(True)
    plt.show()


def HohmannTransfer(r1, r2):
    a = (r1 + r2) / 2  # m :: Major semi-axis of the ellipse in the Hohmann transfer
    t_Hohmann = 1 / 2 * np.sqrt((a ** 3 * 4 * np.pi ** 2) / (G_const * M))  # s :: Time elapsed during the transfer
    ecc = (r2 - r1) / (r2 + r1)  # Eccentricity of the ellipse

    # Inner orbit
    x0, y0, v0x, v0y = Runge_Kutta(f, t_max, dt, w0)


    # Prepare the state for Hohmann speed increase no 1
    v0x[-1] *= np.sqrt(1 + ecc)
    v0y[-1] *= np.sqrt(1 + ecc)
    w1 = np.array([x0[-1], v0x[-1], y0[-1], v0y[-1]])

    print("Inner orbit calculated")

    # Elliptical orbit
    x1, y1, v1x, v1y = Runge_Kutta(f,  t_Hohmann, dt, w1)

    # Prepare the state for Hohmann speed increase no 2
    v1x[-1] *= 1 / np.sqrt(1 - ecc)
    v1y[-1] *= 1 / np.sqrt(1 - ecc)
    w2 = np.array([x1[-1], v1x[-1], y1[-1], v1y[-1]])

    print("Elliptical orbit calculated")

    # Outer orbit
    x2, y2, v2x, v2y = Runge_Kutta(f, t_max * (r2 / r1) ** (3 / 2), dt, w2)

    print("Outer orbit calculated")

    plt.plot(0, 0, 'ro', label=r"Jorden", color = 'blue')
    plt.plot(x0, y0, label=r"Indre bane LEO", lw=1.5, color = 'red')
    plt.plot(x1, y1, label=r"Hohmann-bane", lw=1.5, color = 'yellow')
    plt.plot(x2, y2, label=r"Ytre bane GEO", lw=1.5, color = 'green')
    plt.axis('equal')
    plt.xlabel(r"$x$ [m]", size=24)
    plt.ylabel(r"$y$ [m]", size=24)
    plt.legend(loc=1)
    plt.grid()
    plt.show()

    x = np.concatenate((x0, x1, x2))
    y = np.concatenate((y0, y1, y2))
    vx = np.concatenate((v0x, v1x, v2x))
    vy = np.concatenate((v0y, v1y, v2y))

    # Calculate energies
    V = -G_const * M * m / np.sqrt(x ** 2 + y ** 2)
    K = 1 / 2 * m * (vx ** 2 + vy ** 2)
    E = V + K

    # Plot energy
    timelist = np.linspace(0, len(E) * dt, len(E))

    plt.plot(timelist, E, label=r" E for satelitten i Hohmann-bane", lw=2)
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'$E$ [J]')
    plt.legend(loc=1)
    plt.grid()
    plt.show()

def SpeedIncrease(r1, r2):
    if r1 == r2:
        return 0
    if r2 < r1:
        r1, r2 = r2, r1

    # Eccentricity
    ecc = np.abs((r2 - r1) / (r2 + r1))
    # Store halvakse
    a = (r1 + r2) / 2
    # Tid i hohmann-bane
    t = 1 / 2 * np.sqrt((a ** 3 * 4 * np.pi ** 2) / (G_const * M))
    print(t, 'tid i Hohmann-bane')
    v0 = np.sqrt(G_const * M / r1)
    dv = v0 * np.sqrt(1 + ecc)
    w0 = np.array([r1, 0, 0, dv])
    x, y, vx, vy = Runge_Kutta(f, t, dt, w0)

    v = np.sqrt(vx[-1] ** 2 + vy[-1] ** 2)
    dv += v / np.sqrt(1 - ecc)

    return dv


def SpeedSteps(N):
    v = 0
    radi = np.linspace(r_leo, r_geo, N + 1)
    for i in range(1, N + 1):
        v += SpeedIncrease(radi[i - 1], radi[i])
    return v


def plotHohmannSpeeds(number):
    x = np.arange(1, number + 1)
    y = np.zeros(len(x))
    for i in range(number):
        y[i] = SpeedSteps(x[i])

    # Plot
    plt.scatter(2 * x, y, label=r"Hohmann-bane ",  s=50,
                linewidths=1.5)
    plt.xlabel(r" $\Delta v_i$", size=24)
    plt.ylabel(r" $\Sigma\Delta v_i$ [m/s]", size=24)
    plt.legend(loc=2)
    plt.grid(True)
    plt.show()

#HohmannTransfer(r_leo, r_geo)
#plotHohmannSpeeds(4)
#SpeedIncrease(r_leo,r_geo)
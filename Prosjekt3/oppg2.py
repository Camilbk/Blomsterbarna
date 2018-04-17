import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#Set common figure parameters:
newparams = {
    'figure.figsize': (16, 5), 'axes.grid': True,
    'lines.linewidth': 1.5, 'font.size': 19, 'lines.markersize' : 10,
    'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral'}
plt.rcParams.update(newparams)

AU = 1.5 * 10 ** 11      # Aphelion distance [m]
yr = 3.2 * 10 ** 7       # Orbital period [s]
G_const = 6.673e-11       # Gravitational constant [N (m/kg)^2]
M = 1.9891e30       # Solar mass [kg]
m = 2.4e23       # Mercury mass [kg]
alpha = 1.1e-8
C = (G_const * M * yr ** 2) / (AU ** 3)
#print(C)

t_max = 0.2435549219  # t_max Corresponds to one single complete orbit
dt = 0.0002
N = int(t_max / dt)  # Number of steps




def Runge_Kutta(X0,V0):
    X_4RK = np.zeros(N)
    Y_4RK = np.zeros(N)
    U_4RK = np.zeros(N)
    V_4RK = np.zeros(N)
    E_4RK = np.zeros(N - 1)  # Total energy
    K_4RK = np.zeros(N - 1)  # Kinetic energy
    P_4RK = np.zeros(N - 1)  # Potential energy

    total_velocity = np.zeros(N - 1)

    time_vector = np.linspace(0, t_max, N - 1)  # in order to plot the energy as a function of time

    def F(X_, Y_, U_, V_):  # dX/dt
        return U_

    def G(X_, Y_, U_, V_):  # dY/dt
        return V_

    def H(X_, Y_, U_, V_):  # dU/dt
        return -C * (X_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3))*(1+alpha/(np.sqrt(X_ ** 2 + Y_ ** 2)) ** 2)

    def I(X_, Y_, U_, V_):  # dV/dt
        return -C * (Y_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3))*(1+alpha/(np.sqrt(X_ ** 2 + Y_ ** 2)) ** 2)

    X_4RK[0] = X0
    V_4RK[0] = V0

    for n in range(N - 1):
        k_x1 = dt * F(X_4RK[n], Y_4RK[n], U_4RK[n], V_4RK[n])
        k_y1 = dt * G(X_4RK[n], Y_4RK[n], U_4RK[n], V_4RK[n])
        k_u1 = dt * H(X_4RK[n], Y_4RK[n], U_4RK[n], V_4RK[n])
        k_v1 = dt * I(X_4RK[n], Y_4RK[n], U_4RK[n], V_4RK[n])

        k_x2 = dt * F(X_4RK[n] + k_x1 / 2, Y_4RK[n] + k_y1 / 2, U_4RK[n] + k_u1 / 2, V_4RK[n] + k_v1 / 2)
        k_y2 = dt * G(X_4RK[n] + k_x1 / 2, Y_4RK[n] + k_y1 / 2, U_4RK[n] + k_u1 / 2, V_4RK[n] + k_v1 / 2)
        k_u2 = dt * H(X_4RK[n] + k_x1 / 2, Y_4RK[n] + k_y1 / 2, U_4RK[n] + k_u1 / 2, V_4RK[n] + k_v1 / 2)
        k_v2 = dt * I(X_4RK[n] + k_x1 / 2, Y_4RK[n] + k_y1 / 2, U_4RK[n] + k_u1 / 2, V_4RK[n] + k_v1 / 2)

        k_x3 = dt * F(X_4RK[n] + k_x2 / 2, Y_4RK[n] + k_y2 / 2, U_4RK[n] + k_u2 / 2, V_4RK[n] + k_v2 / 2)
        k_y3 = dt * G(X_4RK[n] + k_x2 / 2, Y_4RK[n] + k_y2 / 2, U_4RK[n] + k_u2 / 2, V_4RK[n] + k_v2 / 2)
        k_u3 = dt * H(X_4RK[n] + k_x2 / 2, Y_4RK[n] + k_y2 / 2, U_4RK[n] + k_u2 / 2, V_4RK[n] + k_v2 / 2)
        k_v3 = dt * I(X_4RK[n] + k_x2 / 2, Y_4RK[n] + k_y2 / 2, U_4RK[n] + k_u2 / 2, V_4RK[n] + k_v2 / 2)

        k_x4 = dt * F(X_4RK[n] + k_x3, Y_4RK[n] + k_y3, U_4RK[n] + k_u3, V_4RK[n] + k_v3)
        k_y4 = dt * G(X_4RK[n] + k_x3, Y_4RK[n] + k_y3, U_4RK[n] + k_u3, V_4RK[n] + k_v3)
        k_u4 = dt * H(X_4RK[n] + k_x3, Y_4RK[n] + k_y3, U_4RK[n] + k_u3, V_4RK[n] + k_v3)
        k_v4 = dt * I(X_4RK[n] + k_x3, Y_4RK[n] + k_y3, U_4RK[n] + k_u3, V_4RK[n] + k_v3)

        X_4RK[n + 1] = X_4RK[n] + k_x1 / 6 + k_x2 / 3 + k_x3 / 3 + k_x4 / 6
        Y_4RK[n + 1] = Y_4RK[n] + k_y1 / 6 + k_y2 / 3 + k_y3 / 3 + k_y4 / 6
        U_4RK[n + 1] = U_4RK[n] + k_u1 / 6 + k_u2 / 3 + k_u3 / 3 + k_u4 / 6
        V_4RK[n + 1] = V_4RK[n] + k_v1 / 6 + k_v2 / 3 + k_v3 / 3 + k_v4 / 6

        E_4RK[n] = 0.5 * (U_4RK[n + 1] ** 2 + V_4RK[n + 1] ** 2) - C / np.sqrt(X_4RK[n + 1] ** 2 + Y_4RK[n + 1] ** 2)
        K_4RK[n] = 0.5 * m * (U_4RK[n + 1] ** 2 + V_4RK[n + 1] ** 2) ** 2
        P_4RK[n] = -(G_const * M * m / np.sqrt(X_4RK[n + 1] ** 2 + Y_4RK[n + 1] ** 2))
        total_velocity[n] = np.sqrt(U_4RK[n + 1] ** 2 + V_4RK[n + 1] ** 2)

    if (t_max == 0.2435549219 ): #t_max = when Mercury has finished one orbit
        # Find offset
        print("\nSmall offset indicates closed orbit")
        print("Offset in 'x': %0.3e - %0.7e = %0.7e" % (X_4RK[0], X_4RK[-1], X_4RK[N - 1] - X_4RK[0]))
        print("Offset in 'y': %0.3e - %0.7e = %0.7e" % (Y_4RK[0], Y_4RK[-1], Y_4RK[N - 1] - Y_4RK[0]))
        print("Total offset: %0.3e" % np.sqrt((X_4RK[0] - X_4RK[-1]) ** 2 + (Y_4RK[0] - Y_4RK[-1]) ** 2))

        deltaTheta = (np.arctan(Y_4RK[-1]/X_4RK[-1])-np.arctan(Y_4RK[0]/X_4RK[0]))*(360/(2*np.pi))
        print(deltaTheta)


        # Find perihelion seperation:
        r_perihelion = abs(min(Y_4RK))
        print("\nThe parahelion seperation is %0.3f, compared to 0.967." % r_perihelion)

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

    return X_4RK,Y_4RK



# Initial conditions
MX0 = 0.47 #[AU]
MV0 = 8.2 #[AU/yr]
Mx,My = Runge_Kutta(MX0,MV0)


EX0 = 1
EV0 = 2*np.pi
Ex,Ey = Runge_Kutta(EX0,EV0)

plt.figure()
plt.title(r'Perihelion of Mercury, $\alpha = %e$' %alpha)
plt.plot(Mx, My, 'g', label = r'Mercury ($\epsilon=0.206$)',color = 'purple', linewidth = 0.5)
plt.plot(Mx[0],My[0], 'ro', color = 'yellow')
plt.plot(Mx[-1],My[-1],'ro', color = 'blue')
plt.plot([0],[0],'ro', label='Sun')
plt.plot(Ex, Ey, 'g', label = r'Earth ($\epsilon=0$)')
plt.legend()
plt.xlabel(r"$x$ [AU]")
plt.ylabel(r"$y$ [AU]")
plt.axes().set_aspect('equal','datalim')
plt.grid()
plt.show()


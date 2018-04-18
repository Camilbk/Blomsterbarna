import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#Set common figure parameters:
newparams = {
    'figure.figsize': (16, 5), 'axes.grid': True,
    'lines.linewidth': 1.5, 'font.size': 19, 'lines.markersize' : 10,
    'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral'}
plt.rcParams.update(newparams)


Er = 6371e3         #Earth radius [m]
Leo = 322e3         #Low earth orbit [m]
Geo = 35680e3       #Geostationary earth orbit [m]
r = Leo+Er          #Distance from center of earth to satellite [m]
p = 120*60          #Orbital period [s]
G_const = 6.673e-11 #Gravitational constant [N (m/kg)^2]
M = 5.9726e24       #Earth mass [kg]
m = 720             #Satellite mass [kg]
C = (G_const * M * p **2 /r**3)
#print(C)

t_max = 1  # t_max Corresponds to one single complete orbit [s]
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
        return -C * (X_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3))

    def I(X_, Y_, U_, V_):  # dV/dt
        return -C * (Y_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3))

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
SX0 = 1 #[r]
SV0 = 2*np.pi #[m/s]
Sx,Sy = Runge_Kutta(SX0,SV0)

plt.figure()
plt.title('LEO satellite orbiting earth')
plt.plot(Sx, Sy, 'g', label = r'LEO satellite ($\epsilon=0.0$)',color = 'purple', linewidth = 0.5)
plt.plot([0],[0],'ro', label='Earth')
plt.legend()
plt.xlabel(r"$x$ [r]")
plt.ylabel(r"$y$ [r]")
plt.axes().set_aspect('equal','datalim')
plt.grid()
plt.show()

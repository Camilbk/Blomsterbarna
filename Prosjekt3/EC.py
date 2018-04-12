import matplotlib.pyplot as plt
import numpy as np

#Simulerer en planetbane ved bruk av Runge-Kutta

AU  = 1.5e11  # Aphelion distance [m]
yr  = 3.2e7   # Orbital period [s]
G_const = 6.673e-11       # Gravitational constant [N (m/kg)^2]
M = 1.9891e30       # Solar mass [kg]
m = 5.9726e24       # Earth mass [kg]

C =  G_const*M*yr**2/AU**3


t_max = 1
dt = 0.0002
N=int(t_max/dt) #Number of steps
time_vector = np.linspace(0,t_max,N-1) #in order to plot the energy as a function of time

R_EC = np.zeros(N)
X_EC = np.zeros(N)
Y_EC = np.zeros(N)
U_EC = np.zeros(N)
V_EC = np.zeros(N)
E_EC = np.zeros(N - 1)  # Total energy
K_EC = np.zeros(N - 1)  # Kinetic energy
P_EC = np.zeros(N - 1)  # Potential energy

total_velocity = np.zeros(N - 1)

# Initial conditions
X0 = 1
U0 = 0
Y0 = 0
V0 = 2 * np.pi


def Euler_Cromer(t):

    X_EC[0] = X0
    V_EC[0] = V0


    for n in range (N-1):
        R_EC[n] = (X_EC[n]**2+Y_EC[n]**2)**(1/2)
        U_EC[n+1] = U_EC[n]-(C*X_EC[n]/R_EC[n])*dt
        V_EC[n+1] = V_EC[n]-(C*Y_EC[n]/R_EC[n])*dt
        X_EC[n+1] = X_EC[n] + U_EC[n+1]*dt
        Y_EC[n+1] = Y_EC[n] + V_EC[n+1]*dt

        E_EC[n] = 0.5 * (U_EC[n + 1] ** 2 + V_EC[n + 1] ** 2) - C / np.sqrt(X_EC[n + 1] ** 2 + Y_EC[n + 1] ** 2)
        K_EC[n] = 0.5 * m * (U_EC[n + 1] ** 2 + V_EC[n + 1] ** 2) ** 2
        P_EC[n] = -(G_const * M * m / np.sqrt(X_EC[n + 1] ** 2 + Y_EC[n + 1] ** 2))
        total_velocity[n] = np.sqrt(U_EC[n + 1] ** 2 + V_EC[n + 1] ** 2)

    total_E = np.sum(E_EC)/(N-1)
    return total_E
'''''
Euler_Cromer(dt)
plt.figure()
plt.title('Euler-Cromer method')
plt.plot(X_EC, Y_EC, 'g', [0], [0], 'ro')
plt.xlabel(r"$x$ [AU]")
plt.ylabel(r"$y$ [AU]")
plt.grid()
plt.show()

plt.figure()
plt.title('Energy as function of time')
plt.plot(time_vector,E_EC)
plt.xlabel("Time [yr]")
plt.ylabel("Energy [J]")
plt.grid()
plt.show()

plt.figure()
plt.title('Potenial energy as a function of time')
plt.plot(time_vector,P_EC)
plt.ylabel("Energy [J]")
plt.xlabel("Time [yr]")
plt.grid()
plt.show()

plt.figure()
plt.title('Kinetic energy as a function of time')
plt.plot(time_vector,K_EC)
plt.ylabel("Energy [J]")
plt.xlabel("Time [yr]")
plt.grid()
plt.show()

plt.figure()
plt.title('Speed as a function of time')
plt.plot(time_vector,total_velocity)
plt.ylabel("Speed [AU/yr]")
plt.xlabel("Time [yr]")
plt.grid()
plt.show()
'''





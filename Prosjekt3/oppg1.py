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
m = 5.9726e24       # Earth mass [kg]

C = (G_const * M * yr ** 2) / (AU ** 3)
#print(C)

t_max = 1  # t_max=1 Corresponds to one single complete orbit
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

    def F(X_, Y_, U_, V_):  # dX/dtau
        return U_

    def G(X_, Y_, U_, V_):  # dY/dtau
        return V_

    def H(X_, Y_, U_, V_):  # dU/dtau
        return -C * (X_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3))

    def I(X_, Y_, U_, V_):  # dV/dtau
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


    return X_4RK,Y_4RK

    # If t_max was set to exactly one period, it will be interesting
    # to investigate whether or not the orbit is closed (planet returns
    # to its starting position)
    '''
    if (t_max == 1):
        # Find offset
        print("\nSmall offset indicates closed orbit")
        print("Offset in 'x': %0.3e - %0.7e = %0.7e" % (X_4RK[0], X_4RK[-1], X_4RK[N - 1] - X_4RK[0]))
        print("Offset in 'y': %0.3e - %0.7e = %0.7e" % (Y_4RK[0], Y_4RK[-1], Y_4RK[N - 1] - Y_4RK[0]))
        print("Total offset: %0.3e" % np.sqrt((X_4RK[0] - X_4RK[-1]) ** 2 + (Y_4RK[0] - Y_4RK[-1]) ** 2))

        # Find perihelion seperation:
        r_perihelion = abs(min(Y_4RK))
        print("\nThe parahelion seperation is %0.3f, compared to 0.967." % r_perihelion)
    '''



def Euler_Cromer(X0,V0):
    R_EC = np.zeros(N)
    X_EC = np.zeros(N)
    Y_EC = np.zeros(N)
    U_EC = np.zeros(N)
    V_EC = np.zeros(N)
    E_EC = np.zeros(N - 1)  # Total energy
    K_EC = np.zeros(N - 1)  # Kinetic energy
    P_EC = np.zeros(N - 1)  # Potential energy

    total_velocity = np.zeros(N - 1)

    time_vector = np.linspace(0, t_max, N - 1)  # in order to plot the energy as a function of time

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

    return 0


def Keplers_third_law():
    #The orbital period squared is equal to a constant times the semimajor axis cubed
    a_p = np.array([0.39,0.72,1.0,1.52,5.20,9.54,19.19,30.06,39.53]) #Semimajor axis for all planets, Mercury to pluto
    T_p = np.zeros(len(a_p))
    for i in range (len(a_p)):
        T_p[i]=np.sqrt(a_p[i]**3)
        print('T',T_p[i],'  ','a_p',a_p)

    #Plotting log(T) against log(a)
    plt.figure()
    plt.plot(np.log(a_p),np.log(T_p),'ro')
    plt.plot(np.log(a_p),np.log(T_p))
    plt.title('Keplers third law for all planets from Mercury to Pluto')
    plt.xlabel('Semimajor axis log(a) [AU]')
    plt.ylabel('Period log(T) [yr]')
    plt.show()


#Keplers_third_law()

# Initial conditions
EX0 = 1
EV0 = 2*np.pi
Ex,Ey = Runge_Kutta(EX0,EV0)

plt.figure()
plt.title('4th order Runge-Kutta')
plt.plot(Ex, Ey, 'g', label = r'Earth ($\epsilon=0$)')
plt.plot([0],[0],'ro', label='Sun')
plt.legend()
plt.xlabel(r"$x$ [AU]")
plt.ylabel(r"$y$ [AU]")
plt.axes().set_aspect('equal','datalim')
plt.grid()
plt.show()





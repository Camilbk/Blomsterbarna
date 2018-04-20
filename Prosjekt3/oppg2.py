import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Set common figure parameters:
newparams = {
    'figure.figsize': (16, 5), 'axes.grid': True,
    'lines.linewidth': 1.5, 'font.size': 19, 'lines.markersize': 10,
    'mathtext.fontset': 'stix', 'font.family': 'STIXGeneral'}
plt.rcParams.update(newparams)

AU = 1.5 * 10 ** 11  # Aphelion distance [m]
yr = 3.2 * 10 ** 7  # Orbital period [s]
G_const = 6.673e-11  # Gravitational constant [N (m/kg)^2]
M = 1.9891e30  # Solar mass [kg]
m = 2.4e23  # Mercury mass [kg]
C = (G_const * M * yr ** 2) / (AU ** 3)

t_max = 0.2435549219 + .1  # t_max Corresponds to one single complete orbit
dt = 0.0002
N = int(t_max / dt)  # Number of steps


def Runge_Kutta(X0, V0, alpha):
    X_4RK = np.zeros(N)
    Y_4RK = np.zeros(N)
    U_4RK = np.zeros(N)
    V_4RK = np.zeros(N)
    E_4RK = np.zeros(N - 1)  # Total energy
    K_4RK = np.zeros(N - 1)  # Kinetic energy
    P_4RK = np.zeros(N - 1)  # Potential energy

    total_velocity = np.zeros(N - 1)
    radius = np.zeros(N)

    time_vector = np.linspace(0, t_max, N - 1)  # in order to plot the energy as a function of tid

    def F(X_, Y_, U_, V_):  # dX/dt
        return U_

    def G(X_, Y_, U_, V_):  # dY/dt
        return V_

    def H(X_, Y_, U_, V_):  # dU/dt
        return -C * (X_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3)) * (1 + alpha / (np.sqrt(X_ ** 2 + Y_ ** 2)) ** 2)

    def I(X_, Y_, U_, V_):  # dV/dt
        return -C * (Y_ / ((np.sqrt(X_ ** 2 + Y_ ** 2)) ** 3)) * (1 + alpha / (np.sqrt(X_ ** 2 + Y_ ** 2)) ** 2)

    X_4RK[0] = X0
    V_4RK[0] = V0
    radius[0] = np.sqrt(V0 ** 2 + 0 ** 2)
    tid = 0
    precession = 0
    trigger = False
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
        radius[n + 1] = np.sqrt(X_4RK[n + 1] ** 2 + Y_4RK[n + 1] ** 2)
        tid = tid + dt

        if trigger == False and radius[n] < radius[n + 1]:
            trigger = True
            print('n =', n)
        if trigger == True and radius[n] > radius[n + 1]:
            print('n =', n)
            precession = np.arctan(Y_4RK[n] / X_4RK[n])
            print('precession pr Mercury orbit =', precession)
            print('Time =', tid, 'yr')
            break


    #Plots
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

    return X_4RK,Y_4RK,precession



# Initial conditions
alpha = 1.1e-8
MX0 = 0.47 #[AU]
MV0 = 8.2 #[AU/yr]
Mx,My,theta= Runge_Kutta(MX0,MV0,alpha)


EX0 = 1
EV0 = 2*np.pi
Ex,Ey,theta = Runge_Kutta(EX0,EV0,alpha)

plt.figure()
plt.title(r'Percession of Mercury, $\alpha = %e$' %alpha)
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



#Record percession as function of alpha
def percession():
    alpha_vector = np.array([1.1e-8,1.1e-7,1.1e-6,1.1e-5,1.1e-4,1.1e-3,1.1e-2])
    print('alphavec',alpha_vector)
    percession_of_Mercury = np.zeros(len(alpha_vector))
    for a in range (len(alpha_vector)): #run for each alpha
        alpha2 = alpha_vector[a]
        values = Runge_Kutta(0.47,8.2,alpha2) #Calculate orbit of Mercury
        percession_of_Mercury[a] = values[2]
        print('alpha = ',alpha2)
        print("\nThe parahelion seperation is %0.3f"% percession_of_Mercury[a],'theta',values[2])


    # plot per as func of alp
    plt.plot(alpha_vector,percession_of_Mercury)
    plt.title('The percession of Mercury as function of alpha')
    plt.xlabel(r'alpha [$AU^2$]')
    plt.ylabel('Percession of Mercury [degrees]')
    plt.xlim(0,1.1e-2)
    plt.show()

percession()

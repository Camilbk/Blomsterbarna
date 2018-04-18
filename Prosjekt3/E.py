import numpy as np
from matplotlib import pyplot as plt

#Constants

AU = 1.5*10**11
yr = 3.2*10**7
M_e = 6*10**24
M_s = 333480*M_e
M_m = 0.1074*M_e
dt = 0.00002
n = int(10/dt)
a_e = 1
e_e = 0.017
a_m = 1.5237
e_m = 0.0934
R_ae = a_e*(1+e_e)

def Euler_Cromer(x1, y1, vx, vy, dt):
    x = x1 + vx*dt
    y = y1 + vy*dt
    return x, y
A = []
t = []
V = []
x = [1]
y = [0]
Vx = 0
Vy = 2*np.pi # np.sqrt(4*np.pi**2)*np.sqrt(((1-e_e)/a_e*(1+e_e))*(1+M_e/M_s))
r = [0]*n
for i in range(n):
    r[i] = (x[i]**2 + y[i]**2)**0.5
    Vx = Vx - (4*np.pi**2*x[i]*dt)/(r[i]**3)
    Vy = Vy - (4*np.pi**2*y[i]*dt)/(r[i]**3)
    V.append(np.sqrt(Vx**2 + Vy**2))
    x_temp, y_temp = Euler_Cromer(x[i], y[i], Vx, Vy, dt)
    x.append(x_temp)
    y.append(y_temp)
    A.append(r[i]*V[i]*dt/2)
    t.append(dt*i)


R_p = min(r)
R_a = max(r)
print('Difference between perihelion and aphelion: ', R_a-R_p)

plt.figure(1)
plt.plot(x, y, linewidth=1)
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.title("Earth's orbit")
plt.grid()
plt.axis('equal')
# plt.xlim(-1.01, 1.01)
# plt.ylim(-1.01, 1.01)

circle1=plt.Circle((0,0),.1,color='gold')
plt.gcf().gca().add_artist(circle1)

plt.figure(2)
plt.plot(t, A)
plt.xlabel('t (years)')
plt.ylabel('A (AU^2)')
plt.title('Area swept by interconnecting line')

plt.show()

# Initial conditions

x_e = [1]
y_e = [0]
Vx_e = 0
Vy_e = np.sqrt(4*np.pi**2)*np.sqrt(((1-e_e)/a_e*(1+e_e))*(1+M_e/M_s))
x_m = [1.5237]
y_m = [0]
Vx_m = 0
Vy_m = np.sqrt(4*np.pi**2)*np.sqrt(((1-e_m)/a_m*(1+e_m))*(1+M_m/M_s))


for i in range(n):
    r_e = (x_e[i]**2 + y_e[i]**2)**0.5
    r_m = (x_m[i] ** 2 + y_m[i] ** 2) ** 0.5
    r_em = ((x_e[i]-x_m[i]) ** 2 + (y_e[i] - y_m[i]) ** 2) ** 0.5

    Vx_e = Vx_e - (4*np.pi**2*x_e[i]*dt)/(r_e**3) - (4*np.pi**2*(M_m/M_s)*x_e[i]*dt)/(r_em**3)
    Vy_e = Vy_e - (4 * np.pi ** 2 * y_e[i] * dt) / (r_e ** 3) - (4 * np.pi ** 2 * (M_m / M_s) * y_e[i] * dt) / (r_em ** 3)

    Vx_m = Vx_m - (4 * np.pi ** 2 * x_m[i] * dt) / (r_m ** 3) - (4 * np.pi ** 2 * (M_e / M_s) * x_m[i] * dt) / (
    r_em ** 3)
    Vy_m = Vy_m - (4 * np.pi ** 2 * y_m[i] * dt) / (r_m ** 3) - (4 * np.pi ** 2 * (M_e / M_s) * y_m[i] * dt) / (
    r_em ** 3)

    x_etemp, y_etemp = Euler_Cromer(x_e[i], y_e[i], Vx_e, Vy_e, dt)
    x_mtemp, y_mtemp = Euler_Cromer(x_m[i], y_m[i], Vx_m, Vy_m, dt)

    x_e.append(x_etemp)
    y_e.append(y_etemp)

    x_m.append(x_mtemp)
    y_m.append(y_mtemp)




plt.figure(3)
plt.plot(x_e, y_e, linewidth=0.2, label="Earth's orbit")
plt.plot(x_m, y_m, linewidth=0.2, label="Mars' orbit")
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.title('Planetary orbits over 10 years')
plt.legend()
plt.grid()
plt.axis('equal')
# plt.xlim(-1.01, 1.01)
# plt.ylim(-1.01, 1.01)

circle1 = plt.Circle((0,0),.1,color='gold')
plt.gcf().gca().add_artist(circle1)

plt.show()
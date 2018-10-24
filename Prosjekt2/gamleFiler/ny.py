import scipy
import numpy as np
import matplotlib.pyplot as plt
import math
import timeit
import random


# from energy import getEnergy
# from energy import U_ij

def makeGrid(N):  # nxn matrise
    grid = np.zeros((N + 2, N + 2)).astype(np.int16)
    n = len(grid)
    grid[int(np.round(n / 2)), int(np.round(n / 2 - N / 2)): int(np.round(n / 2 - N / 2)) + N] = np.linspace(1, N,
                                                                                                             N).astype(
        np.int16)
    return grid


def findX(grid, x):
    pos = np.argwhere(grid == x)[0]
    col = pos[1]
    row = pos[0]
    return row, col


def rigid_rot(grid, x, lengde):
    if (x > math.floor((lengde / 2) + 1)):
        rot = grid.copy()
        rot[rot <= x] = 0

        rigid = grid.copy()
        rigid[rigid > x] = 0
    else:
        rot = grid.copy()
        rot[rot >= x] = 0

        rigid = grid.copy()
        rigid[rigid < x] = 0

    return rot, rigid


def isLegalTwist(twistedgrid, grid):
    if np.count_nonzero(grid) == np.count_nonzero(twistedgrid):
        bol = True
    else:
        bol = False
    # print('bol = ', bol)
    return bol


def twist(grid, lengde, T):
    twist = False
    kb = 1.38 * np.exp(-23)
    B = 1 / (kb * T)

    while (twist == False):
        x = np.random.randint(1, lengde + 1, None)
        print(x, 'x')

        n = np.random.randint(2, size=None)  # genererer clockwise
        print(n, 'n')

        rot, rigid = rigid_rot(grid, x, lengde)

        row_nonzero = np.count_nonzero(np.count_nonzero(rot, axis=1))  # teller rader med tall i
        col_nonzero = np.count_nonzero(np.count_nonzero(rot, axis=0))  # teller kolonner med tall i

        if row_nonzero > col_nonzero:
            side = row_nonzero
        else:
            side = col_nonzero

        row, col = findX(rigid, x)

        twister = rot[(row - side):(row + side + 1), (col - side):(col + side + 1)]
        twister = np.rot90(twister, (2 * n + 1))

        rot[(row - side):(row + side + 1), (col - side):(col + side + 1)] = twister
        twisted_matrix = np.add(rot, rigid)

        twist = isLegalTwist(twisted_matrix, grid)

    E2 = getEnergy(twisted_matrix, U_ij(15))
    E1 = getEnergy(grid, U_ij(15))
    r = np.random.uniform(0, 1)

    if E2 <= E1:
        return twisted_matrix
    elif r < np.exp(-B * (E2 - E1)):
        return twisted_matrix
    else:
        return grid


def twist_execute(antall_twists, lengde, grid, T):
    for i in range(np.round(antall_twists)):
        print('løkken har gjørt ', i, 'ganger')
        temp = grid
        twisted_matrix = twist(temp, lengde, T)
        grid = twisted_matrix

    return grid


################ OPPGAVE 2 ###############

def U_ij(lengde):
    U_ij = np.zeros((lengde + 1, lengde + 1))

    for i in range(1, (lengde + 1)):
        for j in range(1, lengde + 1):
            np.set_printoptions(precision=2)
            U_ij[i, j] = np.random.uniform(low=(-3.47 * 10 ** (-21)), high=(-10.4 * 10 ** (-21)))
            np.fill_diagonal(U_ij, 0)

    for z in range(1, lengde):
        U_ij[z, z] = 0
        U_ij[z + 1, z] = 0
        U_ij[z, z + 1] = 0

    return U_ij


def nearestNeighbours(grid, U, row, col):
    n = len(grid)
    if row + 1 < n:
        E = + U[grid[row + 1, col], grid[row, col]]
    if row - 1 >= 0:
        E = E + U[grid[row - 1, col], grid[row, col]]
    if col + 1 < n:
        E = E + U[grid[row, col + 1], grid[row, col]]
    if col - 1 >= 0:
        E = E + U[grid[row, col - 1], grid[row, col]]
    return E


# dobbel for løkke, ytre løkke der T går fra 0 til 1500 K, og indreløkke kjører d(T) twists ved hver temp.

def getEnergy(polymer, U):  # finner energien til hele polymeret
    totalEnergy = 0
    for x in range(1, 16):  # Går igjennom lengden på polymeret
        row, col = findX(polymer, x)
        NeighbourEnergy = nearestNeighbours(polymer, U, row, col)
        totalEnergy += NeighbourEnergy
    return totalEnergy


def plotEnergy(grid):
    start = timeit.timeit()

    U = U_ij(15)
    E = np.zeros(1500)
    Temp = np.zeros(1500)

    for T in range(1, 1500, 500):
        antall_twists = np.floor(11000 * np.exp(-0.0015 * (T - 0.9999999))).astype(int)
        polymer = twist_execute(antall_twists, 15, grid, T)
        E = np.append(E, (getEnergy(polymer, U)))
        Temp = np.append(Temp, (T - 0.9999999))
        print(T)

    plt.plot(Temp, E)
    plt.title(r'Gjennomsnittsenergi, $\langle E \rangle$, som funksjon av temperatur, $T$')
    plt.xlabel(r'$T$')
    plt.ylabel(r'$\langle E \rangle$')
    plt.legend()
    plt.savefig('plotEnergy.pdf')
    plt.show()

    end = timeit.timeit()
    print('Tid plotEnergy: ', (end - start) / 60, ' minutter\n')


def plotBindingEnergy(grid):
    start = timeit.timeit()

    U = U_ij(15)
    E = np.zeros(5000)
    twists = np.zeros(5000)
    temp = 0

    for twist in range(1, 5000):
        polymer = twist_execute(twist, 15, grid, temp)
        E = np.append(E, (getEnergy(polymer, U)))
        twists = np.append(twists, twist)
        print(twist)


    plt.plot(twists, E)
    plt.title(r'Bindingsenergi, $E$, som funksjon av antall tvister med $T=0$ K') # TITTEL
    plt.xlabel('Antall tvister')
    plt.ylabel(r'$E$')
    plt.legend()
    plt.savefig('plotBindingEnergy_0K.pdf') # NAVN
    plt.show()

    end = timeit.timeit()
    print('\nTid plotBindingEnergy: ', (end - start) / 60, ' minutter\n')





########## OPPGAVE 3 ##########



########## OPPGAVE 4 ##########

def gradualCooling(grid):
    #number_of_twists = np.linspace(0, 30000, 50) #definerer x-aksen

    E = np.zeros(50)
    lengde_polymer = 15
    U = U_ij(lengde_polymer)


    for T in range (1500,-30): #Finner energien for temperaturer som synker med 30K  fra 1500K til 0K
        polymer=twist_execute(600,lengde_polymer,grid,T) #Tvister polymeret 600 ganger for hver temperatur
        E = np.append(E.getEnergy(polymer,U))
        twist = np.append(twist,600)


    plt.plot(twist, E)
    plt.title(r'Energi, $\langle E \rangle$ ved gradvis kjøling av protein, som funksjon av antall tvister $T$')
    plt.xlabel(r'Antall tvister')
    plt.ylabel(r'$\langle E \rangle$')
    plt.legend()
    plt.savefig('gradualCooling.pdf')
    plt.show()




def main():
    polymer = makeGrid(15)
    plotBindingEnergy(polymer)


main()












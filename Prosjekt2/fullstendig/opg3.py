import numpy as np
from oppgave124 import *
import math
import matplotlib.pyplot as plt
import time

# From T = 0 to T = 1500, step 20

def plotMeanDiameter():
    start_t = time.time()
    T = np.linspace(0.0001, 1500, 20)
    meanDiameter15 = np.zeros(20)
    meanDiameter30 = np.zeros(20)
    dmax15 = 15000
    dmax30 = 30000
    teller = 0

    #for hver temp skal det regnes en diameter
    for temp in T:
        teller += 1
        print(temp, 'temp løkke')
        antall_twists15 = int(np.floor(dmax15 * np.exp(-0.0015 * (temp - 0.999999))))
        print(antall_twists15, 'twists for 15')
        antall_twists30 =  int(np.floor(dmax30 * np.exp(-0.0015 * (temp - 0.999999))))
        print(antall_twists30, 'twists for 30 ')
        polymer15, diameter15 = twist_executeDiameter(antall_twists15, 15, makeGrid(15), temp )  #lager ny polymer for hver temperatur
        polymer30, diameter30 = twist_executeDiameter(antall_twists30, 30, makeGrid(30), temp)  # lager ny polymer for hver temperatur
        meanDiameter15[teller -1] = diameter15
        meanDiameter30[teller -1] = diameter30

    end_t = time.time()
    plt.figure(1)
    plt.plot(T, meanDiameter15, lw=2, label=r"15 monomers", color="b")
    plt.plot(T, meanDiameter30, lw=2, label=r"30 monomers", color="crimson")
    plt.xlabel(r"$T$  [K]")
    plt.ylabel(r"$\langle L \rangle$ ")
    plt.legend(loc="best")
    plt.grid()
    plt.figure(2)
    plt.plot(T, meanDiameter30, lw=2, label=r"30 monomers", color="crimson")
    plt.xlabel(r"$T$  [K]")
    plt.ylabel(r"$\langle L \rangle$ ")
    plt.legend(loc="best")
    plt.grid()
    plt.figure(3)
    plt.plot(T, meanDiameter15, lw=2, label=r"15 monomers", color="b")
    plt.xlabel(r"$T$  [K]")
    plt.ylabel(r"$\langle L \rangle$ ")
    plt.legend(loc="best")
    plt.grid()
    print((start_t-end_t)/60 ,'min')
    plt.show()



def calculateDiameter(grid, lengde):
    longest = 0
    for i in range(1, lengde):
        row1, col1 = findX(grid, i)
        for j in range(i, lengde):
            row2, col2 = findX(grid, j)
            current = math.sqrt((row2 - row1) ** 2 + (col2 - col1) ** 2)
            if current > longest:
                longest = current
    return longest


def twist_executeDiameter(antall_twists,lengde,grid, T):
    meanDiameter = np.zeros(np.round(antall_twists))
    for i in range(np.round(antall_twists)):
        temp = grid
        twisted_matrix = twist(temp, lengde, T)
        grid = twisted_matrix
        meanDiameter[i] = calculateDiameter(grid, lengde)
    diameter = np.average(meanDiameter)
    return grid, diameter


def twist(grid, lengde,T):
    twist = False
    kb = 1.38*10**(-23)
    B = 1/(kb * T)


    while(twist == False):
        x = np.random.randint(1,lengde + 1,None)
        n = np.random.randint(2, size = None) #genererer clockwise

        rot, rigid = rigid_rot(grid, x, lengde)

        row_nonzero = np.count_nonzero(np.count_nonzero(rot, axis=1)) #teller rader med tall i
        col_nonzero = np.count_nonzero(np.count_nonzero(rot, axis=0)) #teller kolonner med tall i

        if row_nonzero > col_nonzero:
            side = row_nonzero
        else:
            side = col_nonzero

        row, col = findX(rigid, x)


        twister = rot[(row-side):(row+side+1),(col-side):(col+side+1)]
        twister = np.rot90(twister,(2*n+1))

        rot[(row - side):(row + side + 1), (col - side):(col + side + 1)] = twister
        twisted_matrix = np.add(rot,rigid)

        twist = isLegalTwist(twisted_matrix,grid)

    E2 = getEnergy(twisted_matrix, U_ij(lengde),lengde )
    E1 = getEnergy(grid, U_ij(lengde), lengde)
    r = np.random.uniform(0,1)

    if E2 <= E1:
        return twisted_matrix
    elif r < np.exp(-B*(E2-E1)):
        return twisted_matrix
    else:
        return grid


def getEnergy(polymer, U, lengde): #finner energien til hele polymeret
    totalEnergy = 0
    for x in range (1, lengde+1): #Går igjennom lengden på polymeret
        row, col = findX(polymer,x )
        NeighbourEnergy = nearestNeighbours(polymer,U,row,col)
        totalEnergy += NeighbourEnergy
    return totalEnergy*0.5

def nearestNeighbours(grid, U, row, col ):
    n = len(grid)
    if row + 1 < n:
        E =+ U[grid[row + 1, col], grid[row, col]]
    if row - 1 >= 0:
        E =+ U[grid[row - 1, col], grid[row, col]]
    if col + 1 < n:
        E =+ U[grid[row, col + 1], grid[row, col]]
    if col -1 >= 0:
        E =+ U[grid[row, col - 1], grid[row, col]]
    return E

def U_ij(lengde):
    U_ij = np.zeros((lengde+1, lengde+1))

    for i in range (1,(lengde+1)):
         for j in range (1,lengde +1):
            np.set_printoptions(precision=2)
            U_ij[i,j] = np.random.uniform(low = (-3.47*10**(-21)), high= (-10.4*10**(-21)))
            np.fill_diagonal(U_ij,0)

    for z in range (1,lengde):
            U_ij[z,z] = 0
            U_ij[z+1,z] = 0
            U_ij[z,z+1] = 0

    return U_ij


def isLegalTwist(twistedgrid, grid):
    if np.count_nonzero(grid) == np.count_nonzero(twistedgrid):
        bol = True
    else:
        bol = False
    #print('bol = ', bol)
    return bol


plotMeanDiameter()

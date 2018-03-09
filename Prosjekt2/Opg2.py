import scipy
import numpy as np
import matplotlib.pyplot as plt
import math
import timeit

from ny import twist_execute
from ny import findX
from ny import makeGrid

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

def nearestNeighbours(grid, U, row, col ):
    n = len(grid)
    if row + 1 < n:
        E =+ U[grid[row + 1, col], grid[row, col]]
    if row - 1 >= 0:
        E = E + U[grid[row - 1, col], grid[row, col]]
    if col + 1 < n:
        E = E + U[grid[row, col + 1], grid[row, col]]
    if col -1 >= 0:
        E = E + U[grid[row, col - 1], grid[row, col]]
    return E
# dobbel for løkke, ytre løkke der T går fra 0 til 1500 K, og indreløkke kjører d(T) twists ved hver temp.

def getEnergy(polymer, U): #finner energien til hele polymeret
    totalEnergy = 0
    for x in range (1,16): #Går igjennom lengden på polymeret
        row, col = findX(polymer,x )
        NeighbourEnergy = nearestNeighbours(polymer,U,row,col)
        totalEnergy += NeighbourEnergy
    return totalEnergy

def plotEnergy(grid):
    start = timeit.timeit()

    U = U_ij(15)
    E = np.zeros(1500)
    Temp = np.zeros(1500)

    for T in range (1,1500,1):
        antall_twists = np.floor(11000 * np.exp(-0.0015 * (T-0.9999999))).astype(int)
        polymer = twist_execute(antall_twists,15, grid)
        E = np.append(E,(getEnergy(polymer,U)))
        Temp = np.append(Temp, (T-0.9999999))
        print(T)


    plt.plot(Temp, E)
    plt.show()

    end = timeit.timeit()
    print(end-start)

def plotBindingEnergy(grid):
    start = timeit.timeit()

    U = U_ij(15)
    E = np.zeros(5000)
    Twists = np.zeros(5000)

    for Twist in range(1, 5000):
        polymer = twist_execute(Twist, 15, grid)
        E = np.append(E, (getEnergy(polymer, U)))
        Twists = np.append(Twists, Twist)
        print(Twist)

    plt.plot(Twists, E)
    plt.show()

    end = timeit.timeit()
    print(end - start)



def main():
    polymer = makeGrid(15)
    plotEnergy(polymer)
    plotBindingEnergy(polymer)

    '''#antall_twists = np.floor(11000 * np.exp(-0.0015 * 10)).astype(int)
    antall_twists = 5
    grid = twist_execute(antall_twists, 15)
    print(grid)
    U = U_ij(15)
    print(U)
    energy = getEnergy(grid,U)
    print(energy)'''



main()

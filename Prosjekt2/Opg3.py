import numpy as np
from Vitber.Prosjekt2.ny import *
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
        antall_twists30 = int(np.floor(dmax30 * np.exp(-0.0015 * (temp - 0.999999))))
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

plotMeanDiameter()

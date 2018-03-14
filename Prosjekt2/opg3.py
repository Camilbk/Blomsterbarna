import numpy as np
from Vitber.Prosjekt2.ny import *
import matplotlib.pyplot as plt

# From T = 0 to T = 1500, step 20

def plotMeanDiameter():
    start = timeit.timeit()
    T = np.linspace(0, 1500, 20)
    meanDiameter15 = np.zeros(20)
    meanDiameter30 = np.zeros(20)
    dmax15 = 11000
    dmax30 = 21000

    #for hver temp skal det regnes en diameter
    for temp in T:
        antall_twists15 = int(np.floor(dmax15 * np.exp(-0.0015 * (T-0.999999))))
        antall_twists30 = int(np.floor(dmax30 * np.exp(-0.0015 * (T - 0.999999))))
        polymer15, diameter15 = twist_executeDiameter(antall_twists15, 15, makeGrid(15), temp )  #lager ny polymer for hver temperatur
        polymer30, diameter30 = twist_executeDiameter(antall_twists30, 30, makeGrid(30), temp)  # lager ny polymer for hver temperatur
        meanDiameter15[temp] = diameter15
        meanDiameter30[temp] = diameter30

    plt.plot(T, meanDiameter15, lw=3, label=r"15 monomers", color="crimson")
    plt.plot(T, meanDiameter30, lw=3, label=r"30 monomers", color="crimson")
    plt.xlabel(r"$T$  [K]", size=20)
    plt.ylabel(r"$\langle L \rangle$  [1]", size=20)
    plt.legend(loc="best")
    plt.grid()
    plt.show()



def calculateDiameter(grid):


def twist_executeDiameter(antall_twists,lengde,grid, T):
    meanDiameter = np.zeros(np.round(antall_twists))
    for i in range(np.round(antall_twists)):
        temp = grid
        twisted_matrix = twist(temp, lengde, T)
        grid = twisted_matrix
        meanDiameter[i] = calculateDiameter(grid)
    diameter = np.average(meanDiameter)
    return grid, diameter


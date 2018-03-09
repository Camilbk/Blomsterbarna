from Vitber.Prosjekt2.Opg1 import findX
import scipy
import numpy as np
import matplotlib.pyplot as plt
import math

from ny import makeGrid

'''
Define a straight polymer consisting of 15 monomers. Follow the example in Sec. 1.4, and plot the mean energy, ⟨E⟩,
as function of temperature. Choose the relevant parameters as in Sec. 1.3. That is, choose Uij randomly
between −3.47· 10−21 and −10.4· 10−21 (in units of J), and vary T between 0 K and 1500 K. You need to do d(T)
twists at each temperature T, where d(T) = dmax exp{−sT}, as explained in Sec. 1.4. Choose dmax and s so that you
obtain a sufficiently convergent5 result,.
'''


# Lag polymermatrise

'''
Denote the potential be- tween the amino acids i and j at positions ri and rj as U(ri − rj), 
and the expression for the mean energy would then follow as ⟨E⟩ = sum (states (protein structure)) = sum U(ri − rj)
 That is, if A1 is nearest neighbour to A4, a weak bond forms between these two monomers
'''

# Lag nabomatrise U, med tall tilfeldig på intervall  (−3.47· 10−21 ,   −10.4· 10−21) (in units of J). (we define Ui j to be zero if |i − j| ≤ 1)
# Polymer og nabo må kobles sammen slik at hvis twist så må nabo endres.

def getEnergy(polymer, U): #finner energien til hele polymeret
    totalEnergy = 0
    for x in range (1,16): #Går igjennom lengden på polymeret
        row, col = findX(polymer,x )
        NeighbourEnergy = nearestNeighbours()
        totalEnergy += NeighbourEnergy
    return totalEnergy


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


def U_ij(lengde):
    U_ij = np.zeros((lengde+1, lengde+1))
    print(U_ij)

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


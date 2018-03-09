import scipy
import numpy as np
import matplotlib.pyplot as plt
import math

from energy import getEnergy
from energy import U_ij

def makeGrid(N):  # nxn matrise
    grid = np.zeros((N+2, N+2)).astype(np.int16)
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

def rigid_rot(grid,x,lengde):
    if (x > np.floor(lengde/2)+1):
        rot = grid.copy()
        rot[rot <= x]=0

        rigid=grid.copy()
        rigid[rigid > x]=0
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
    #print('bol = ', bol)
    return bol

def twist(grid, lengde,T):
    twist = False
    kb = 1.38 * np.exp(-23)
    B = 1 / (kb * T)

    while(twist == False):
        rot, rigid = rigid_rot(grid, x, lengde)

        x = np.random.randint(1, high=lengde + 1)
        n = np.random.randint(2, size=1)  # clockwise

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

    E2 = getEnergy(twisted_matrix, U_ij(15))
    E1 = getEnergy(grid, U_ij(15))
    r = np.random.uniform(0,1)

    if E2 <= E1:
        return twisted_matrix
    elif r < np.exp(-B*(E2-E1)):
        return twisted_matrix
    else:
        return grid

def twist_execute(antall_twists,lengde,twisted_matrix, T):
    for i in range(np.round(antall_twists)):
        grid = twisted_matrix
        polymer = twist(grid, lengde, T)
        twisted_matrix = polymer

    return twisted_matrix




def main():
    polymer = makeGrid(15)
    twist(polymer,15,500)




















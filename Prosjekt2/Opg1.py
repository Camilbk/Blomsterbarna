import scipy
import numpy as np
import math


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
        bol = 1
    else:
        bol = 0
    #print('bol = ', bol)
    return bol

def twist(grid,x,lengde):
    isLegalTwist = False
    # x = np.randint #SKAL FINNE RANDOM TALL AV POLYMER

    while(not isLegalTwist):
        rot,rigid=rigid_rot(grid,x,lengde)

        n = np.random.randint(2, size=1)  # clockwise

        row_nonzero = np.count_nonzero(np.count_nonzero(rot, axis=1)) #teller rader med tall i
        col_nonzero = np.count_nonzero(np.count_nonzero(rot, axis=0)) #teller kolonner med tall i

        if row_nonzero > col_nonzero:
            side = row_nonzero
        else:
            side = col_nonzero

        row,col =findX(rigid,x)


        twister = rot[(row-side):(row+side+1),(col-side):(col+side+1)]
        twister = np.rot90(twister,(2*n+1))

        rot[(row - side):(row + side + 1), (col - side):(col + side + 1)] = twister
        twisted_matrix = np.add(rot,rigid)

        isLegalTwist = isLegalTwist(twisted_matrix,grid)


    return twisted_matrix



def main():
    lengde=11
    grid=makeGrid(lengde)
    #grid[grid > 5] = 0
    #print(grid)
    print(twist(grid,4,11))




    #rot=grid.copy()
   # print(rot)
   # rot[rot > 4] = 0
   # print(rot)

main()
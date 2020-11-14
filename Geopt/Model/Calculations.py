import math
from ase import Atom
import numpy as np


def ExtraSpace(atomsList, periods, xmlList):
    # Some atoms are much larger than others so the cell must be adjusted for this.
    size = 0
    lastAtom = ''
    for eachAtom in atomsList:
        # Don't look up the same element multiple times
        if eachAtom != lastAtom:
            for i in range(periods):
                for j in xmlList[i]:
                    if j[0].text == eachAtom:
                        group = int(j[4].text)
                        relativeSize = (18 - group) * i
                        # Update size to find largest atom
                        if relativeSize > size:
                            size = relativeSize
                            print(size)
        lastAtom = eachAtom
    return size/10


def EvenSpacing(atomsList, extraSpace):
    # Space out atoms evenly in a 3D box.
    numItems = len(atomsList)
    # Find the cube root because a box is a cube. Round up with ceiling. This will be the dimensions of the box.
    # Add 2 units so that we can centre the atoms easily later on.
    axis = math.ceil(numItems ** (1 / 3))
    axisSquared = axis ** 2
    dimensions = axis + 2 + extraSpace*axis
    boxSize = (dimensions, dimensions, dimensions)
    print("boxSize:", boxSize)
    print("segments per plane:", axis)
    print("extraSpace:", extraSpace)
    # the box is split into segments of size 1,1,1 plus extra spacing for the largest atom
    coordinates = []
    for i in range(axis):
        for j in range(axis):
            for k in range(axis):
                # Use axis as the numerical base to get the index of the atom list without another loop.
                index = i*axisSquared + j*axis + k
                # Create a matrix for the co-ordinates of each atom in the molecule
                coordinates.append(Atom(atomsList[index], [i+i*extraSpace, j+j*extraSpace, k+k*extraSpace]))
                # Don't create too many co-ordinates. numAtoms is unlikely to divide perfectly into the grid.
                if len(coordinates) == numItems:
                    return boxSize, coordinates








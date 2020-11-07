import math


def EvenSpacing(atomsList):
    # Space out atoms evenly in a 3D box.
    numItems = len(atomsList)
    # Find the cube root because a box is a cube. Round up with ceiling. This will be the dimensions of the box.
    # Add 2 units so that we can centre the atoms easily later on.
    axis = math.ceil(numItems ** (1 / 3))
    axisSquared = axis ** 2
    dimensions = axis + 2
    boxSize = (dimensions, dimensions, dimensions)
    # the box is split into segments of size 1,1,1
    coordinates = []
    for i in range(axis):
        for j in range(axis):
            for k in range(axis):
                # Use axis as the numerical base to get the index of the atom list without another loop.
                index = i*axisSquared + j*axis + k
                print(index)
                coordinates.append((atomsList[index], [i, j, k]))
                # Don't create too many co-ordinates. numAtoms is unlikely to divide perfectly into the grid.
                if len(coordinates) == numItems:
                    print(coordinates)
                    print(len(coordinates))
                    return boxSize, coordinates


testList = ('C'*18)
EvenSpacing(testList)

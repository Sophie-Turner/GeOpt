import math
from ase import Atom
from ase.data import atomic_numbers, covalent_radii
from Model.InteractWithData import GetXML


def ExtraSpace(atomsList):
    # Some atoms are much larger than others so the cell must be adjusted for this.
    treerootMain, treerootF = GetXML()
    maxSize = 0
    covRads = []
    covRad = None
    lastAtom = ''
    for eachAtom in atomsList:
        # Don't look up the same element multiple times.
        if eachAtom != lastAtom:
            # Most elements are in the main block so check there first.
            isFound, maxSize, covRad = FindAtom(eachAtom, treerootMain, 7, maxSize)
            if isFound is False:
                # Check in the F block if it hasn't been found.
                isFound, maxSize, covRad = FindAtom(eachAtom, treerootF, 2, maxSize)
        lastAtom = eachAtom
        covRads.append(covRad)
    return maxSize/16, covRads


def FindAtom(atomToFind, xmlList, periods, maxSize):
    isFound = False
    covRad = None
    for i in range(periods):
        for j in xmlList[i]:
            if j[0].text == atomToFind:
                isFound = True
                group = int(j[4].text)
                relativeSize = (18 - group) * i
                covRad = covalent_radii[atomic_numbers[atomToFind]]
                # Update size to find largest atom
                if relativeSize > maxSize:
                    maxSize = relativeSize
    return isFound, maxSize, covRad


def EvenSpacing(atomsList, extraSpace):
    # Space out atoms evenly in a 3D box.
    numItems = len(atomsList)
    # Find the cube root because a box is a cube. Round up with ceiling. This will be the dimensions of the box.
    # Add 2 units so that we can centre the atoms easily later on.
    axis = math.ceil(numItems ** (1 / 3))
    axisSquared = axis ** 2
    # The box is split into segments of size 1,1,1 plus extra spacing for the largest atom.
    dimensions = axis + 2 + extraSpace*axis
    boxSize = (dimensions, dimensions, dimensions)
    # The coordinates list is kept separate so it can be altered easily during the EA.
    coordinates = []
    # This list is for creating the initial model and it's quicker to build it up during this loop.
    atomObjectList = []
    for i in range(axis):
        x = i + i * extraSpace
        for j in range(axis):
            y = j + j * extraSpace
            for k in range(axis):
                # Use axis as the numerical base to get the index of the atom list without another loop.
                index = i*axisSquared + j*axis + k
                # Create a matrix for the co-ordinates of each atom in the molecule
                z = k+k*extraSpace
                coordinates.append([x, y, z])
                atomObjectList.append(Atom(atomsList[index], [x, y, z]))
                # Don't create too many co-ordinates. numAtoms is unlikely to divide perfectly into the grid.
                if len(coordinates) == numItems:
                    return boxSize, coordinates, atomObjectList

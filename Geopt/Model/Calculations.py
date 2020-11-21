import math
from ase import Atom
import numpy as np
from Model.InteractWithData import GetXML


def ExtraSpace(atomsList):
    # Some atoms are much larger than others so the cell must be adjusted for this.
    treerootMain, treerootF = GetXML()
    maxSize = 0
    lastAtom = ''
    for eachAtom in atomsList:
        # Don't look up the same element multiple times.
        if eachAtom != lastAtom:
            # Most elements are in the main block so check there first.
            isFound, maxSize = FindAtom(eachAtom, treerootMain, 7, maxSize)
            if isFound is False:
                # Check in the F block if it hasn't been found.
                isFound, maxSize = FindAtom(eachAtom, treerootF, 2, maxSize)
        lastAtom = eachAtom
    return maxSize/10


def FindAtom(atomToFind, xmlList, periods, maxSize):
    isFound = False
    for i in range(periods):
        for j in xmlList[i]:
            if j[0].text == atomToFind:
                isFound = True
                print("The atom was found.")
                group = int(j[4].text)
                relativeSize = (18 - group) * i
                # Update size to find largest atom
                if relativeSize > maxSize:
                    maxSize = relativeSize
    return isFound, maxSize


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
    print("boxSize:", boxSize)
    print("segments per plane:", axis)
    print("extraSpace:", extraSpace)
    # The coordinates list is kept separate so it can be altered easily during the EA.
    coordinates = []
    # This list of for creating the initial model and it's quicker to build it up during this loop.
    atomObjectList = []
    for i in range(axis):
        for j in range(axis):
            for k in range(axis):
                # Use axis as the numerical base to get the index of the atom list without another loop.
                index = i*axisSquared + j*axis + k
                # Create a matrix for the co-ordinates of each atom in the molecule
                x = i+i*extraSpace
                y = j+j*extraSpace
                z = k+k*extraSpace
                coordinates.append([x, y, z])
                atomObjectList.append(Atom(atomsList[index], [x, y, z]))
                # Don't create too many co-ordinates. numAtoms is unlikely to divide perfectly into the grid.
                if len(coordinates) == numItems:
                    return boxSize, coordinates, atomObjectList


def SortByEnergy(population):
    # Sort the population by energy values.
    print("length of population: ", len(population))
    indexSmall = 0
    indexBig = len(population) - 1

    print("Before the sort: ", population)



def Traverse(indexSmall, indexBig, population):
    pivot = round(indexBig / 2)
    pivotValue = population[pivot][2]
    print("pivot index: ", pivot)
    print("pivot value: ", pivotValue)
    while indexSmall < indexBig:
        while population[indexSmall][2] < pivotValue:
            indexSmall = indexSmall + 1
            print("indexSmall incremented to ", indexSmall)
        while population[indexBig][2] > pivotValue:
            indexBig = indexBig - 1
            print("indexBig decremented to ", indexBig)
        smallMolecule = population[indexSmall]
        bigMolecule = population[indexBig]
        if smallMolecule[2] > pivotValue > bigMolecule[2]:
            Swap(indexSmall, indexBig, population)
        elif smallMolecule[2] > pivotValue:
            Swap(indexSmall, pivot, population)
        elif pivotValue > bigMolecule[2]:
            Swap(pivot, indexBig, population)
    print("After this iteration: ", population)
    


def Swap(smallIndex, bigIndex, population):
    print(population[smallIndex][2], "is bigger than", population[bigIndex][2], ". Swapping.")
    tempCopy = population[bigIndex]
    population[bigIndex] = population[smallIndex]
    population[smallIndex] = tempCopy






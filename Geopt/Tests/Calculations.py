import math
from Model.InteractWithData import GetXML


def ExtraSpace(atomsList, periods, xmlList):
    # Some atoms are much larger than others so the cell must be adjusted for this.
    size = 0
    lastAtom = 'C'
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


def EvenSpacing(atomsList):
    # Space out atoms evenly in a 3D box.
    numItems = len(atomsList)
    # Find the cube root because a box is a cube. Round up with ceiling. This will be the dimensions of the box.
    # Add 2 units so that we can centre the atoms easily later on.
    axis = math.ceil(numItems ** (1 / 3))
    axisSquared = axis ** 2
    dimensions = axis + 2
    # the box is split into segments of size 1,1,1
    boxSize = (dimensions, dimensions, dimensions)
    print("boxSize:", boxSize)
    # The coordinates list is kept separate so it can be altered easily during the EA.
    coordinates = []
    # This list of for creating the initial model and it's quicker to build it up during this loop.
    atomObjectList = []
    for i in range(axis):
        for j in range(axis):
            for k in range(axis):
                # Use axis as the numerical base to get the index of the atom list without another loop.
                index = i*axisSquared + j*axis + k
                print(index)
                coordinates.append([i, j, k])
                atomObjectList.append((atomsList[index], [i, j, k]))
                # Don't create too many co-ordinates. numAtoms is unlikely to divide perfectly into the grid.
                if len(coordinates) == numItems:
                    print(coordinates)
                    print(len(coordinates))
                    return boxSize, coordinates, atomObjectList


def SortByEnergy(population):
    # Sort the population by energy values.
    print("length of population: ", len(population))
    indexSmall = 0
    indexBig = len(population) - 1
    print("Before the sort: ", population)
    Traverse(indexSmall, indexBig, population)


def Traverse(indexSmall, indexBig, population):
    pivot = round((indexBig - indexSmall) / 2)
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


testList = ('C', 'H', 'Li', 'Na')
EvenSpacing(testList)
treerootMain, treerootF = GetXML()
ExtraSpace(testList, 7, treerootMain)

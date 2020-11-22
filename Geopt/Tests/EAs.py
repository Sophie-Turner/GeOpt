from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model.Populations import Population
from numpy import random


def StartEA(elementsList):
    # Set up and initialise our template molecule to start with.
    thisPopulation = Population(elementsList)
    boxSize = thisPopulation.boxSize
    firstCoordinates = thisPopulation.initPositions
    atomObjectList = thisPopulation.initAtomsObject
    numAtoms = len(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()
    parentEnergy = GetEnergy(parentMolecule)
    print("Parent energy: ", parentEnergy)
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/parent.png", parentMolecule,
          rotation='10x,30y,0z')

    if numAtoms == 1:
        # If there's only 1 atom we can skip all this...
        pass

    # Start with some permutations of an initial model.
    population = [[parentMolecule, firstCoordinates, parentEnergy]]
    for i in range(numAtoms * 2):
        MakeNewMolecule(elementsList, firstCoordinates, 1000, None, boxSize, population, False, False, True)

    Evolve(elementsList, boxSize, firstCoordinates, population)


def Evolve(elementsList, boxSize, firstCoordinates, population):
    # Ranges of random atom movements.
    width = boxSize[0]
    changeSizes = [width/30, width/20, width/16, width/12, width/8, width/4]
    # See how many iterations it takes.
    iterations = 0
    similarity = 0
    lastBestEnergy = 1000
    # End if the best energy doesn't change much for several consecutive iterations.
    while similarity < 10:
        lastBestEnergy = population[0][2]
        # Selection.
        population = RankByE(population, 2)
        bestCoordinates = population[0][1]
        # New child molecules.
        for i in range(len(elementsList) * 2):
            # Make some permutations.
            MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, None, boxSize, population, False, False, True)
            # Make some with random mutations.
            MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, changeSizes[0], boxSize, population, 2, False, False)
            # Introduce some random strangers.
            MakeNewMolecule(elementsList, firstCoordinates, lastBestEnergy, None, boxSize, population, 3, False, False)
            # Try some crossovers
            MakeNewMolecule(elementsList, None, lastBestEnergy, changeSizes[0], boxSize, population, 2, True, False)
            # Move Hydrogens.
            MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, None, boxSize, population, 5, False, False)


        # Update the stopping criterion.
        newBestEnergy = population[0][2]
        if abs(lastBestEnergy - newBestEnergy) < 0.01:
            similarity = similarity + 1
        else:
            similarity = 0
        iterations = iterations + 1

    print("Iterations performed: ", iterations)
    print("The best energy found was: ", lastBestEnergy)
    print("Population size: ", len(population))
    print("Final co-ordinates: ", bestCoordinates)

    # write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/crossover.png", population[65][0],
    #       rotation='10x,30y,0z')
    # write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/mutation.png", population[25][0],
    #        rotation='10x,30y,0z')
    # write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/random.png", population[45][0],
    #        rotation='10x,30y,0z')
    # write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/permutation.png", population[5][0],
    #        rotation='10x,30y,0z')
    # write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/dad.png", population[1][0],
    #       rotation='10x,30y,0z')
    # write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/Hydrogens.png", population[75][0],
    #        rotation='10x,30y,0z')

    RankByE(population, 1)
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/optimised.png", population[0][0],
          rotation='10x,30y,0z')


def MakeNewMolecule(elementsList, inCoordinates, bestE, changeSize, boxSize, population, mutate, cross, permute):
    if cross is True:
        inCoordinates = Crossover(population)
    if permute is True:
        inCoordinates = random.permutation(inCoordinates)
    childCoordinates, atomsObject = PrepareChild(elementsList, inCoordinates)
    if mutate == 1:
        MoveAtomsUniform(atomsObject, changeSize)
    elif mutate == 2:
        MoveAtomsGauss(atomsObject, 0, changeSize)
    elif mutate == 3:
        FillCellGauss(atomsObject, boxSize)
    elif mutate == 4:
        FillCellUniform(atomsObject, boxSize)
    elif mutate == 5:
        MoveHydrogen(atomsObject, boxSize)
    childMolecule = GenerateChild(atomsObject, boxSize)
    childEnergy = GetEnergy(childMolecule)
    # Don't keep any extremely terrible structures.
    if childEnergy >= bestE * 100:
        MakeNewMolecule(elementsList, inCoordinates, bestE, changeSize, boxSize, population, mutate, cross, permute)
    population.append([childMolecule, childCoordinates, childEnergy])


def PrepareChild(elementsList, startCoordinates):
    childAtomsObject = []
    # Need to copy the parent's co-ordinates so we don't alter the parent molecule at the same time.
    childCoordinates = startCoordinates[:]
    for i in range(len(elementsList)):
        childAtomsObject.append(Atom(elementsList[i], childCoordinates[i]))
    return childCoordinates, childAtomsObject


def GenerateChild(childAtomsObject, boxSize):
    child = Atoms(childAtomsObject, cell=boxSize)
    # Translate atoms to the centre of the unit cell. It's OK if the atoms stick out of their cell.
    child.center()
    return child


def GetEnergy(molecule):
    # Set up the ase force calculator for finding energies.
    molecule.calc = EMT()
    energy = molecule.get_potential_energy()
    return energy


def MoveAtomsUniform(itsAtoms, changeSize):
    for eachAtom in itsAtoms:
        eachAtom.position += ((random.randint(-changeSize, changeSize, 3)) / 100)



def MoveAtomsGauss(itsAtoms, mean, sigma):
    for eachAtom in itsAtoms:
        eachAtom.position += (random.normal(mean, sigma, 3))


def FillCellUniform(itsAtoms, boxSize):
    for eachAtom in itsAtoms:
        eachAtom.position = (random.random(3)*boxSize/2)


def FillCellGauss(itsAtoms, boxSize):
    centre = boxSize[0]/2
    for eachAtom in itsAtoms:
        eachAtom.position = (random.normal(centre, centre/3, 3))


def MoveHydrogen(itsAtoms, boxSize):
    # Hydrogens tend to get lost so push them all onto the Carbon atoms.
    change = boxSize[0]/12
    moveTo = []
    carbons = 0
    for eachAtom in itsAtoms:
        if eachAtom.symbol == 'C':
            moveTo.append(eachAtom.position)
        elif eachAtom.symbol == 'H':
            if len(moveTo) > 0:
                try:
                    eachAtom.position = moveTo[carbons] + (random.normal(0, change, 3))
                    carbons += 1
                except:
                    carbons = 0
                    eachAtom.position = moveTo[carbons] + (random.normal(0, change, 3))
            else:
                eachAtom.position += (random.normal(0, change * 2, 3))



def RankByE(population, numToKeep):
    # I originally planned to use a quicksort but decided that since the list was small
    # there was no need to make it complicated.
    rankedPopulation = []
    bestMolecule = None
    while len(rankedPopulation) < numToKeep:
        bestEnergy = 1000
        for eachMember in population:
            if eachMember[2] < bestEnergy:
                bestEnergy = eachMember[2]
                bestMolecule = eachMember
        rankedPopulation.append(bestMolecule)
        population.remove(bestMolecule)
    return rankedPopulation


def Crossover(population):
    mumCoordinates = population[0][1]
    howMany = len(mumCoordinates)
    mumOrDad = random.randint(0, 2, howMany)
    childCoordinates = []
    for i in range(howMany):
        whichParent = mumOrDad[i]
        thisPoint = population[whichParent][1][i]
        childCoordinates.append(thisPoint)
    return childCoordinates







testAcetone = ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H']
testWater = ['H', 'H', 'O']
testMethane = ['C', 'H', 'H', 'H', 'H']
testC12 = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C']
testN2O4 = ['N', 'N', 'O', 'O', 'O', 'O']
testBenzene = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
testAcetonitrile = ['C', 'H', 'H', 'H', 'C', 'N']
testCO2 = ['C', 'C', 'O']
StartEA(testWater)

# How to remove old objects from memory:
#    for eachAtom in child1AtomsObject:
#        print("removing ", eachAtom)
#        del(eachAtom)
#    print("removing ", child1Molecule)
#    del(child1Molecule)

class Optimiser:

    # User can adjust parameters

    def __init__(self, molecule):
        self.molecule = molecule

    def EAtest(self):
        self.molecule.GetEnergy()

from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model.Populations import Population
from numpy import random

# Ranges of random atom movements.
global changeSizes
changeSizes = [10, 20, 50, 100, 150]

def StartEA(elementsList):

    # Set up and initialise our template molecule to start with.
    thisPopulation = Population(elementsList)
    boxSize = thisPopulation.boxSize
    firstCoordinates = thisPopulation.initPositions
    atomObjectList = thisPopulation.initAtomsObject

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()
    parentEnergy = GetEnergy(parentMolecule)
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/parent.png", parentMolecule,
          rotation='10x,30y,0z')

    if len(elementsList) == 1:
        # If there's only 1 atom we can skip all this...
        pass

    # Create 3 permutations from this initial parent.
    population = [[parentMolecule, firstCoordinates, parentEnergy]]
    for i in range(3):
        MakeNewMolecule(elementsList, firstCoordinates, None, boxSize, population, False, False, True)

    Evolve(elementsList, boxSize, firstCoordinates, population)


def Evolve(elementsList, boxSize, firstCoordinates, population):
    # See how many iterations it takes.
    iterations = 0
    similarity = 0
    lastBestEnergy = 1000
    # End if the best energy doesn't change much for 3 consecutive iterations.
    while similarity < 10:
        lastBestEnergy = population[0][2]
        # Selection.
        population = RankByE(population, 1)
        bestCoordinates = population[0][1]
        # New child molecules.
        for i in range(5):
            # Make 3 with random mutations but no crossover.
            MakeNewMolecule(elementsList, bestCoordinates, changeSizes[1], boxSize, population, True, False, False)
            # Make 3 permutations with mutations.
            MakeNewMolecule(elementsList, bestCoordinates, changeSizes[1], boxSize, population, True, False, True)
            # Introduce 2 random strangers.
            MakeNewMolecule(elementsList, firstCoordinates, changeSizes[3], boxSize, population, True, False, True)

        # Update the stopping criterion.
        newBestEnergy = population[0][2]
        if abs(lastBestEnergy - newBestEnergy) < 0.01:
            similarity = similarity + 1
        else:
            similarity = 0
        iterations = iterations + 1
        print("Iteration: ", iterations)
        print("Energy: ", newBestEnergy, "\n")

    print("Iterations performed: ", iterations)
    print("The best energy found was: ", lastBestEnergy)
    print("Population size: ", len(population))
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/optimised.png", population[0][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child1.png", population[1][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child2.png", population[2][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child3.png", population[3][0],
          rotation='10x,30y,0z')


def MakeNewMolecule(elementsList, inCoordinates, changeSize, boxSize, population, mutate, cross, permute):
    if cross is True:
        inCoordinates = Crossover(population)
    if permute is True:
        inCoordinates = random.permutation(inCoordinates)
    childCoordinates, atomsObject = PrepareChild(elementsList, inCoordinates)
    if mutate is True:
        MoveAllAtoms(atomsObject, changeSize)
    childMolecule = GenerateChild(atomsObject, boxSize)
    childEnergy = GetEnergy(childMolecule)
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


def MoveAllAtoms(itsAtoms, changeSize):
    for eachAtom in itsAtoms:
        eachAtom.position += ((random.randint(-changeSize, changeSize, 3)) / 100)


def RankByE(population, numToKeep):
    # I originally planned to use a quicksort but decided that since the list was small
    # there was no need to make it complicated.
    rankedPopulation = []
    bestMolecule = None
    while len(rankedPopulation) < numToKeep:
        bestEnergy = 1000
        for eachMember in population:
            print("This molecule's energy: ", eachMember[2])
            if eachMember[2] < bestEnergy:
                bestEnergy = eachMember[2]
                bestMolecule = eachMember
        rankedPopulation.append(bestMolecule)
        population.remove(bestMolecule)
    print("Ranked population: \n", rankedPopulation)
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


testList = ['C', 'C', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H']
testList2 = ['H', 'H', 'O']
testList3 = ['C', 'H', 'H', 'H', 'H']
StartEA(testList)

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

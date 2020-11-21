from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model import Molecules, Calculations
from numpy import random

def StartEA(elementsList):
    changeSizes = [10, 50, 150]  # Ranges of random atom movements.

    # Set up and initialise our template molecule to start with.
    boxSize, firstCoordinates, atomObjectList = Molecules.SetUpMolecule(elementsList)
    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()
    parentEnergy = GetEnergy(parentMolecule)

    if len(elementsList) == 1:
        # If there's only 1 atom we can skip all this...
        pass

    bestEnergy = parentEnergy
    bestMolecule = parentMolecule

    # Create 3 permutations from this initial parent.
    population = [[parentMolecule, firstCoordinates, parentEnergy]]
    for i in range(3):
        permutedCoordinates = random.permutation(firstCoordinates)
        childCoordinates, childAtomsObject = PrepareChild(elementsList, permutedCoordinates)
        childMolecule = GenerateChild(childAtomsObject, boxSize)
        childEnergy = GetEnergy(childMolecule)
        population.append([childMolecule, childCoordinates, childEnergy])

    # Sort the molecules into order of lowest energy and keep the two best ones.
    population = RankByE(population, 2)
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/dad.png", population[1][0],
          rotation='10x,30y,0z')

    # Make 2 child molecules from these parents
    for i in range(2):
        # Do a random crossover and a small random mutation
        childCoordinates = Crossover(population)
        childCoordinates, childAtomsObject = PrepareChild(elementsList, childCoordinates)
        MoveAllAtoms(childAtomsObject, changeSizes[0])
        childMolecule = GenerateChild(childAtomsObject, boxSize)
        childEnergy = GetEnergy(childMolecule)
        population.append([childMolecule, childCoordinates, childEnergy])

    # Get rid of the dad and introduce a random stranger
    population.pop(1)
    permutedCoordinates = random.permutation(firstCoordinates)
    randCoordinates, randAtomsObject = PrepareChild(elementsList, permutedCoordinates)
    MoveAllAtoms(randAtomsObject, changeSizes[2])
    randMolecule = GenerateChild(randAtomsObject, boxSize)
    randEnergy = GetEnergy(randMolecule)
    population.append([randMolecule, randCoordinates, randEnergy])


    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/mum.png", population[0][0],
         rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child1.png", population[1][0],
         rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child2.png", population[2][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/random.png", population[3][0],
          rotation='10x,30y,0z')


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
    # I originally planned to use a quicksort but decided that since the list was small there was no need to make it complicated.
    rankedPopulation = []
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







testList = ['C', 'O', 'H', 'H', 'H']
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

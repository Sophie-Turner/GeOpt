from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from Model.Populations import Population
from numpy import random
from concurrent import futures


def SetUpVasp():
    calc = Vasp(
        xc='pbe',  # exchange-correlation functional
        nbands=6,  # number of bands
        encut=350,  # planewave cutoff
        ismear=1,  # Methfessel-Paxton smearing
        sigma=0.01,  # very small smearing factor for a molecule
    )
    return calc

def MakeNewMolecule(elementsList, inCoordinates, bestE, changeSize, boxSize, population, mutate, cross, permute, calc):
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
    childEnergy = GetEnergy(childMolecule, calc)
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


def GetEnergy(molecule, calc):
    # Set up the ase force calculator for finding energies.
    molecule.calc = calc
    energy = molecule.get_potential_energy()
    return energy


def MoveOneAtomGauss(fixedAtom, atomToMove, sigma):
    mean = fixedAtom.position
    atomToMove.position = (random.normal(mean, sigma, 3))


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


def PermuteAtoms(itsAtoms, inCoordinates):
    newCoordinates = random.permutation(inCoordinates)
    for eachAtom in itsAtoms:
        eachAtom.position = newCoordinates[eachAtom]


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
testCO2 = ['C', 'O', 'O']
testAlox = ['Al', 'Al', 'O', 'O', 'O']
testEthane = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']
testH2 = ['H', 'H']


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

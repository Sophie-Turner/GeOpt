from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model import Molecules
from numpy import random

def StartEA(elementsList):
    changeSizes = [10, 50, 150]  # Ranges of random atom movements.

    # Set up and initialise our template molecule to start with.
    boxSize, coordinates, atomObjectList = Molecules.SetUpMolecule(elementsList)
    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()
    parentEnergy = GetEnergy(parentMolecule)
    bestEnergy = parentEnergy
    bestMolecule = parentMolecule
    # Create 3 children from this initial parent using large random ranges for mutation.
    childrenList = []
    for i in range(3):
        permutedCoordinates = random.permutation(coordinates)
        childCoordinates, childAtomsObject = PrepareChild(elementsList, permutedCoordinates)
        #MoveAllAtoms(childAtomsObject, changeSizes[2])
        childMolecule = GenerateChild(childAtomsObject, boxSize)
        childEnergy = GetEnergy(childMolecule)
        print("Energy:", childEnergy)
        if childEnergy < bestEnergy:
            bestEnergy = childEnergy
            bestMolecule = childMolecule
        childrenList.append([childMolecule, childCoordinates, childAtomsObject, childEnergy])

    print("Parent energy: ", parentEnergy)
    print("Parent molecule: ", parentMolecule)
    print("Child 1 molecule: ", childrenList[0][0])
    print("Parent atoms object: ", atomObjectList)
    print("Child 1 atoms object: ", childrenList[0][2])
    print("The best energy is", bestEnergy)

    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/parent.png", parentMolecule,
         rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child1.png", childrenList[0][0],
         rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child2.png", childrenList[1][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child3.png", childrenList[2][0],
          rotation='10x,30y,0z')


def PrepareChild(elementsList, parentCoordinates):
    childAtomsObject = []
    childCoordinates = parentCoordinates[:]
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


testList = ['C', 'H', 'H', 'H', 'O']
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

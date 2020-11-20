from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model import Molecules
from numpy import random

def StartEA(elementsList):
    # Set up and initialise our template molecule to start with.
    boxSize, coordinates, atomObjectList = Molecules.SetUpMolecule(elementsList)
    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()

    # Create 3 children from this initial parent using large random ranges for mutation.
    childrenList = []
    for i in range(3):
        childMolecule, childCoordinates, childAtomsObject = GenerateChild(elementsList, coordinates, boxSize)
        childrenList.append([childMolecule, childCoordinates, childAtomsObject])

    parentEnergy = GetEnergy(parentMolecule)
    child1Energy = GetEnergy(childrenList[0][0])

    print("Parent molecule: ", parentMolecule)
    print("Child 1 molecule: ", childrenList[0][0])
    print("Parent atoms object: ", atomObjectList)
    print("Child 1 atoms object: ", childrenList[0][2])
    print("Parent energy: ", parentEnergy)
    print("Child energy: ", child1Energy)

    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/parent.png", parentMolecule,
         rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child1.png", childrenList[0][0],
         rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child2.png", childrenList[1][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child3.png", childrenList[2][0],
          rotation='10x,30y,0z')



def GenerateChild(elementsList, parentCoordinates, boxSize):
    changeSizes = [10, 50, 150]  # Ranges of random atom movements.
    childAtomsObject = []
    childCoordinates = parentCoordinates[:]
    for i in range(len(elementsList)):
        childAtomsObject.append(Atom(elementsList[i], childCoordinates[i]))
    MoveAllAtoms(childAtomsObject, changeSizes[2])
    child = Atoms(childAtomsObject, cell=boxSize)
    # Translate atoms to the centre of the unit cell. It's OK if the atoms stick out of their cell.
    child.center()
    return child, childCoordinates, childAtomsObject



def GetEnergy(molecule):
    # Set up the ase force calculator for finding energies.
    molecule.calc = EMT()
    energy = molecule.get_potential_energy()
    return energy


def MoveAllAtoms(itsAtoms, changeSize):
    for eachAtom in itsAtoms:
        eachAtom.position += ((random.randint(-changeSize, changeSize, 3)) / 100)


testList = ['C', 'H', 'H', 'H', 'H']
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

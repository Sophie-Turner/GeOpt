from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model import Molecules
from numpy import random

def StartEA(elementsList):
    boxSize, coordinates, atomObjectList = Molecules.SetUpMolecule(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentEnergy = GetEnergy(parentMolecule)

    GenerateChild(elementsList, coordinates, boxSize)
    child1Molecule, child1Coordinates, child1AtomsObject = GenerateChild(elementsList, coordinates, boxSize)
    MoveAllAtoms(child1AtomsObject)

    print("Parent molecule: ", parentMolecule)
    print("Child molecule: ", child1Molecule)
    print("Parent atoms object: ", atomObjectList)
    print("Child atoms object: ", child1AtomsObject)

    # How to remove old objects from memory:
    for eachAtom in child1AtomsObject:
        print("removing ", eachAtom)
        del(eachAtom)
    print("removing ", child1Molecule)
    del(child1Molecule)



    #childEnergy = GetEnergy(childMolecule1)

    #print("parent energy: ", parentEnergy, "\nchild energy: ", childEnergy)
    #print("parent positions: \n", parentMolecule.get_positions(), "\nchild positions: \n", childMolecule1.get_positions())
    #write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/parent.png", parentMolecule,
    #      rotation='10x,30y,0z')
    #write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child.png", childMolecule1,
    #      rotation='10x,30y,0z')



def GenerateChild(elementsList, parentCoordinates, boxSize):
    childAtomsObject = []
    childCoordinates = parentCoordinates[:]
    for i in range(len(elementsList)):
        childAtomsObject.append(Atom(elementsList[i], childCoordinates[i]))
    child = Atoms(childAtomsObject, cell=boxSize)
    child.calc = EMT()
    return child, childCoordinates, childAtomsObject



def GetEnergy(molecule):
    # Translate atoms to the centre of the unit cell.
    molecule.center()
    # Set up the ase force calculator for finding energies.
    molecule.calc = EMT()
    energy = molecule.get_potential_energy()
    return energy


def MoveAllAtoms(itsAtoms):
    for eachAtom in itsAtoms:
        eachAtom.position += ((random.randint(-50, 50, 3)) / 100)


testList = ['H', 'H', 'O']
StartEA(testList)


class Optimiser:

    # User can adjust parameters

    def __init__(self, molecule):
        self.molecule = molecule

    def EAtest(self):
        self.molecule.GetEnergy()

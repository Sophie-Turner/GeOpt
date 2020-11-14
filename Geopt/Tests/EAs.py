from ase import Atoms, Atom
from ase.io import write
from ase.calculators.emt import EMT
from Model import Molecules
from numpy import random

def StartEA(elementsList):
    boxSize, atomObjectList = Molecules.SetUpMolecule(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentEnergy = GetEnergy(parentMolecule)
    print("parent atoms: ", atomObjectList)

    MoveAllAtoms(atomObjectList)

    childMolecule = Atoms(atomObjectList, cell=boxSize)
    childEnergy = GetEnergy(childMolecule)
    print("child atoms: ", atomObjectList)

    print("parent energy: ", parentEnergy, "\nchild energy: ", childEnergy)
    print("parent positions: \n", parentMolecule.get_positions(), "\nchild positions: \n", childMolecule.get_positions())
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/parent.png", parentMolecule,
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/child.png", childMolecule,
          rotation='10x,30y,0z')


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

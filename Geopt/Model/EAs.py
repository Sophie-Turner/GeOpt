import numpy as np
from numpy import random
from Model import Molecules
from ase import Atom, Atoms
from ase.calculators.emt import EMT
from ase.io import write
# If there is only 1 atom we do nothing.
# If there are 2 atoms we need to only move one of them and only on one plane.


def StartEA(elementsList):
    boxSize, coordinates, atomObjectList = Molecules.SetUpMolecule(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentEnergy = GetEnergy(parentMolecule)

    child1Molecule, child1Coordinates, child1AtomsObject = GenerateChild(elementsList, coordinates, boxSize)
    childEnergy = GetEnergy(child1Molecule)


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
    # energy in electron volts
    energy = molecule.get_potential_energy()
    return energy


def MoveAllAtoms(itsAtoms):
    for eachAtom in itsAtoms:
        eachAtom.position += ((random.randint(-50, 50, 3)) / 100)




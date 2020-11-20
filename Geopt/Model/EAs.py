import numpy as np
from numpy import random
from Model import Molecules
from ase import Atom, Atoms
from ase.calculators.emt import EMT
from ase.io import write
# If there is only 1 atom we do nothing.
# If there are 2 atoms we need to only move one of them and only on one plane.


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
        if childEnergy < bestEnergy:
            bestEnergy = childEnergy
            bestMolecule = childMolecule
        childrenList.append([childMolecule, childCoordinates, childAtomsObject, childEnergy])


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

import numpy as np
from numpy import random
from Model import Molecules
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import write
# If there is only 1 atom we do nothing.
# If there are 2 atoms we need to only move one of them and only on one plane.


def StartEA(elementsList):
    boxSize, atomObjectList = Molecules.SetUpMolecule(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentEnergy = GetEnergy(parentMolecule)

    MoveAllAtoms(atomObjectList)

    childMolecule = Atoms(atomObjectList, cell=boxSize)
    childEnergy = GetEnergy(childMolecule)


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




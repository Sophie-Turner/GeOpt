import numpy as np
from numpy import random
from Model.Molecules import Molecule
# If there is only 1 atom we do nothing.
# If there are 2 atoms we need to only move one of them and only on one plane.


def StartEA(elementsList):
    # Create an instance of the molecule class and create the molecule.
    thisMolecule = Molecule(elementsList, None, None, None, 'EMT')
    startingPositions = thisMolecule.ModelMolecule()
    oldEnergy = thisMolecule.GetEnergy()
    # Put this inside a loop... and run (4) parallel processes.
    MoveAllAtoms(startingPositions)
    newEnergy = thisMolecule.GetEnergy()
    #print("old energy: ", oldEnergy, "   new energy: ", newEnergy)


def MoveAllAtoms(itsAtoms):
    for eachAtom in itsAtoms:
        MutatePositions(eachAtom)


def MoveOneAtom(itsAtoms, atomToMove):
    pass


def MutatePositions(atomToMove):
    print("old atom", atomToMove)
    atomToMove.position += ((random.randint(-50, 50, 3)) / 100)
    print("new atom", atomToMove)
    return atomToMove




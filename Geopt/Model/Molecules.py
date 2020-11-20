from Model.Calculations import *


def SetUpMolecule(atoms):
    # Find the largest atom's relative size and adjust the cell for this.
    extraSpace = ExtraSpace(atoms)
    # Calculate cell size and positions for atom placement.
    boxSize, atomObjectList = EvenSpacing(atoms, extraSpace)
    return boxSize, atomObjectList












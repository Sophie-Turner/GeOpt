from Model.Calculations import *
from Model.InteractWithData import GetXML


def SetUpMolecule(atoms):
    # Find the largest atom's relative size and adjust the cell for this.
    boxSize = None
    atomObjectList = None
    treerootMain, treerootF = GetXML()
    isValid, extraSpace = ExtraSpace(atoms, 7, treerootMain)
    if isValid is True:
        # Calculate cell size and positions for atom placement.
        boxSize, atomObjectList = EvenSpacing(atoms, extraSpace)
    return boxSize, atomObjectList










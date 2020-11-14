from Model.Calculations import *
from Model.InteractWithData import GetXML


def SetUpMolecule(atoms):
    # Find the largest atom's relative size and adjust the cell for this.
    treerootMain, treerootF = GetXML()
    extraSpace = ExtraSpace(atoms, 7, treerootMain)
    # Calculate cell size and positions for atom placement.
    boxSize, atomObjectList = EvenSpacing(atoms, extraSpace)
    return boxSize, atomObjectList










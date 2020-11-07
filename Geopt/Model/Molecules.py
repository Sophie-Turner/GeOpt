# This is the Atomic Simulation Environment (ASE) package
from ase import Atoms
from ase.io import write
from Model.Calculations import *
from Model.InteractWithData import GetXML

class Molecule:

    def __init__(self, atoms, functionalGroups, ions, charge):
        # The properties a molecule has
        self.atoms = atoms
        self.functionalGroups = functionalGroups
        self.ions = ions
        self.charge = charge
        if self.charge is None:
            self.charge = 0

        self.structure = None


    def ModelMolecule(self):
        # Create a model of this molecule using ASE.
        # Find the largest atom's relative size and adjust the cell for this.
        treerootMain, treerootF = GetXML()
        extraSpace = ExtraSpace(self.atoms, 7, treerootMain)
        # Calculate cell size and positions for atom placement.
        boxSize, atomObjectList = EvenSpacing(self.atoms, extraSpace)
        self.structure = Atoms(atomObjectList, cell=boxSize)
        # Translate atoms to the centre of the unit cell.
        self.structure.center()
        # See if it worked!
        write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/testcell.png", self.structure)




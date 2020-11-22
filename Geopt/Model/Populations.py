from Model.Calculations import *


class Population:
    # Initialise common parameters for a population of molecules.
    def __init__(self, atoms):
        # Find the largest atom's relative size and adjust the cell for this.
        self.extraSpace = ExtraSpace(atoms)
        # Calculate cell size and positions for atom placement.
        self.boxSize, self.initPositions, self.initAtomsObject = EvenSpacing(atoms, self.extraSpace)

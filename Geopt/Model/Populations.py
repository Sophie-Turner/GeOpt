from Model.Calculations import *


class Population:
    # Initialise common parameters for a population of molecules.
    def __init__(self, atoms, mutSize):
        # Find the largest atom's relative size and adjust the cell for this.
        self.extraSpace, self.covRads = ExtraSpace(atoms)
        # Calculate cell size and positions for atom placement.
        boxDims, self.initPositions, self.initAtomsObject = EvenSpacing(atoms, self.extraSpace)
        self.boxSize = []
        for edge in boxDims:
            edge = edge + mutSize
            self.boxSize.append(edge)
        del boxDims

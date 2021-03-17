from Model.Calculations import *


class Population:
    # Initialise common parameters for a population of molecules.
    def __init__(self, atoms, mutSize):
        # Find the largest atom's relative size and adjust the cell for this.
        self.boxSize, self.covRads = [], []
        self.extraSpace, covRads = ExtraSpace(atoms)
        for covRad in covRads:
            covRad = covRad + (mutSize/10)
            self.covRads.append(covRad)
        # Calculate cell size and positions for atom placement.
        boxDims, self.initPositions, self.initAtomsObject = EvenSpacing(atoms, self.extraSpace)
        for edge in boxDims:
            edge = edge + mutSize
            self.boxSize.append(edge)
        del boxDims, covRads

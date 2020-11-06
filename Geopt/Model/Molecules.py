# This is the Atomic Simulation Environment (ASE) package
from ase import Atoms, Atom
from ase.io import write

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
        # Create a model of this molecule using ASE
        atomObjectList = []
        i = 0
        for eachAtom in self.atoms:
            atomObjectList.append(Atom(eachAtom, [i, 0, 0]))
            i=i+0.25
        self.structure = Atoms(atomObjectList, cell=(6, 6, 6))

        # See if it worked!
        self.structure.center()  # translate atoms to centre of unit cell
        write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/testcell.png", self.structure)




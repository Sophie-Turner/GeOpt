class Molecule:

    def __init__(self, atoms, functionalGroups, ions, charge):
        # The properties a molecule has
        self.atoms = atoms
        self.functionalGroups = functionalGroups
        self.ions = ions
        self.charge = charge
        if self.charge is None:
            self.charge = 0


class Isomer:

    # An isomer IS A molecule - class structure...
    def __init__(self):
        pass
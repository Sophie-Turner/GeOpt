from Controller.Shared import *
from tkinter import messagebox

# Long strings that are shown. Stored them here to avoid cluttering up the code.
txtPbc = ('Periodic boundary conditions are used to approximate\n'
          'a larger system from repeating the unit cell.')
txtManyMolecule = ('The many-molecule algorithm is an evolutionary algorithm\n'
                   'which creates a population of versions of the molecule\n'
                   'with different configurations, applies mutations to\n'
                   'generations of versions and selects the best versions\n'
                   'based on the lowest total potential energy. It is\n'
                   'recommended for molecules with three to twelve atoms.')
txtPerAtom = ('The per-atom exhaustive test is an algorithm which moves\n'
              'each atom in the molecule around each other atom,\n'
              'constantly testing for the configuration with the lowest\n'
              'total potential energy. The amount of processing required\n'
              'is proportional to the number of atoms factorial and it\n'
              'is recommended for molecules with two or three atoms.')
txtPopSize = ('The population size is the number of versions of the\n'
              'molecule that will exist at the same time at each step in\n'
              'the algorithm. A larger population makes the algorithm take\n'
              'longer to complete but may offer more potential configurations.')


def ShowTxtPbc(window):
    messagebox.showinfo(parent=window, title='Periodic boundary conditions', message=txtPbc)


def ShowTxtManyMolecule(window):
    messagebox.showinfo(parent=window, title='Many-molecule Evolution', message=txtManyMolecule)


def ShowTxtPerAtom(window):
    messagebox.showinfo(parent=window, title='Per-atom Exhaustive Test', message=txtPerAtom)


def ShowTxtPopSize(window):
    messagebox.showinfo(parent=window, title='Population size', message=txtPopSize)


def EnableEA():
    pass


def DisableEA():
    pass






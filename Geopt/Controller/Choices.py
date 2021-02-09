from Controller.Shared import *
from tkinter import messagebox
from os import cpu_count

# Long strings that are shown. Stored them here to avoid cluttering up the code.
titlePbc = 'Periodic boundary conditions'
txtPbc = ('Periodic boundary conditions are used to approximate\n'
          'a larger system from repeating the unit cell.')
titleManyMolecule = 'Many-molecule Evolution'
txtManyMolecule = ('The many-molecule algorithm is an evolutionary algorithm\n'
                   'which creates a population of versions of the molecule\n'
                   'with different configurations, applies mutations to\n'
                   'generations of versions and selects the best versions\n'
                   'based on the lowest total potential energy. It is\n'
                   'recommended for molecules with three to twelve atoms.')
titlePerAtom = 'Per-atom Exhaustive Test'
txtPerAtom = ('The per-atom exhaustive test is an algorithm which moves\n'
              'each atom in the molecule around each other atom,\n'
              'constantly testing for the configuration with the lowest\n'
              'total potential energy. The amount of processing required\n'
              'is proportional to the number of atoms factorial and it\n'
              'is recommended for molecules with two or three atoms.')
titlePopSize = 'Population size'
txtPopSize = ('The population size is the number of versions of the\n'
              'molecule that will exist at the same time at each step in\n'
              'the algorithm, in each process. A larger population makes\n'
              'the algorithm take longer to complete but may offer more\n'
              'potential configurations.')
titleCores = 'Parallel processes'
txtCores = ('The number of parallel processes is how many times the entire '
            'algorithm will run. It is recommended that this number is equal '
            'to the number of CPU cores in your computer. This is detected '
            'and set as default. A higher number may offer more potential '
            'configurations, but you are strongly advised to avoid setting '
            'this to a higher number than the number of CPU cores for molecules '
            'with more than two atoms.')
titleCoreErr = 'Unable to detect CPU cores'
txtCoreErr = ('Unable to detect the number of CPU cores in your computer.\n'
              'Setting default processes to one.')

messages = [(titleManyMolecule, txtManyMolecule), (titlePerAtom, txtPerAtom), (titlePbc, txtPbc),
            (titlePopSize, txtPopSize), (titleCores, txtCores), (titleCoreErr, txtCoreErr)]


def ShowMessage(window, which):
    message = messages[which]
    messagebox.showinfo(parent=window, title=message[0], message=message[1])






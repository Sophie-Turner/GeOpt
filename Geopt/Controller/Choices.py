from Controller.Shared import *
from tkinter import messagebox
from os import cpu_count

# Long strings that are shown. Stored them here to avoid cluttering up the code.
titlePbc = 'Periodic boundary conditions'
txtPbc = ('Periodic boundary conditions are used to approximate\n'
          'a larger system from repeating the unit cell.')
titleManyMolecule = 'Many-molecule Evolution'
txtManyMolecule = ('The many-molecule algorithm is an evolutionary algorithm '
                   'which creates a population of versions of the molecule '
                   'with different configurations, applies mutations to '
                   'generations of versions and selects the best versions '
                   'based on the lowest total potential energy.\nIt is '
                   'recommended for molecules with three to twelve atoms.')
titlePerAtom = 'Per-atom Exhaustive Test'
txtPerAtom = ('The per-atom exhaustive test is an algorithm which moves '
              'each atom in the molecule around each other atom, '
              'constantly testing for the configuration with the lowest '
              'total potential energy.\nThe amount of processing required '
              'is proportional to the number of atoms factorial and it '
              'is recommended for molecules with two to four atoms.')
titlePopSize = 'Population size'
txtPopSize = ('The population size is the number of versions of the '
              'molecule that will exist at the same time at each step in '
              'the algorithm, in each process.\nA larger population makes '
              'the algorithm take longer to complete but may offer more '
              'potential configurations.')
titleCores = 'Parallel processes'
txtCores = ('The number of parallel processes is how many times the entire '
            'algorithm will run.\nIt is recommended that this number is equal '
            'to the number of CPU cores in your computer. This is detected '
            'and set by default.\nA higher number may offer more potential '
            'configurations, but you are advised to avoid setting this to a '
            'higher number than the number of CPU cores for molecules with '
            'more than two atoms.')
titlePlotLimit = 'Plot data size limit'
txtPlotLimit = ('This is the maximum number of data points in the datasets used to '
                'create plots.\nA smaller number will speed up execution time and '
                'may also create clearer graphs.')
titleCoreErr = 'Unable to detect CPU cores'
txtCoreErr = ('Unable to detect the number of CPU cores in your computer.\n'
              'Setting default processes to one.')

messages = [(titleManyMolecule, txtManyMolecule), (titlePerAtom, txtPerAtom), (titlePbc, txtPbc),
            (titlePopSize, txtPopSize), (titleCores, txtCores), (titlePlotLimit, txtPlotLimit),
            (titleCoreErr, txtCoreErr)]


def ShowMessage(window, which):
    message = messages[which]
    messagebox.showinfo(parent=window, title=message[0], message=message[1])






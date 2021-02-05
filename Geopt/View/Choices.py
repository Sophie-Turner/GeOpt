from tkinter import messagebox
from math import exp
from Controller.Shared import *
from View.Shared import SetUpWindow
from View.Analysis import StartAnalysis

global formula


def ChooseFeatures(elementsList, boxText):
    global formula
    formula = boxText
    window = tk.Toplevel()
    window.attributes('-topmost', 'true')
    SetUpWindow(window)
    # Show the user their chosen formula.
    lblHelp = Label(window, text=boxText, font=('Agency FB', 18), fg='white', bg='#222222')
    lblHelp.pack(pady=10)
    optionsFrame = tk.Frame(window)

    # The place where helpful info will be shown.
    Label(optionsFrame, text='The many-molecule algorithm is an evolutionary algorithm\n'
                             'which creates a population of versions of the molecule\n'
                             'with different configurations, applies mutations to\n'
                             'generations of versions and selects the best versions\n'
                             'based on the lowest total potential energy. It is\n'
                             'recommended for molecules with three to twelve atoms.')\
        .grid(row=0, column=0, rowspan=3)
    Label(optionsFrame, text='The per-atom exhaustive test is an algorithm which moves\n'
                             'each atom in the molecule around each other atom,\n'
                             'constantly testing for the configuration with the lowest\n'
                             'total potential energy. The amount of processing required\n'
                             'is proportional to the number of atoms factorial and it\n'
                             'is recommended for molecules with two or three atoms.')\
        .grid(row=3, column=0, rowspan=3)
    Label(optionsFrame, text='Periodic boundary conditions are used to approximate\n'
                             'a larger system from repeating the unit cell.') \
        .grid(row=6, column=0, rowspan=3)


    # The place where options will be chosen.
    Label(optionsFrame, text='Optimisation algorithm', font=('Agency FB', 14)).grid(row=0, column=1)
    # Set up radio buttons.
    # Algorithms for the user to choose from.
    algos = [('Many-molecule evolution', 0),
             ('Per-atom exhaustive test', 1)]
    algo = tk.IntVar()
    # Default MMEA algo.
    algo.set(0)
    for name, value in algos:
        tk.Radiobutton(optionsFrame, text=name, variable=algo, value=value).grid(row=value+1, column=1)
    # Periodic boundary conditions.
    Label(optionsFrame, text='Periodic boundary conditions', font=('Agency FB', 14)).grid(row=3, column=1)
    pbcs = [('Off', False), ('On', True)]
    pbc = tk.BooleanVar()
    pbc.set(False)
    for name, value in pbcs:
        tk.Radiobutton(optionsFrame, text=name, variable=pbc, value=value).grid(row=value+4, column=1)


    btnBuild = Button(optionsFrame, text='Build molecule', command=(lambda: ProceedToAlgo(elementsList, algo.get(),
                                                                                          pbc.get())),
                      font=('Agency FB', 14))
    btnBuild.grid(row=5, column=2)
    #SetColours(optionsFrame)
    optionsFrame.pack()
    window.mainloop()


def ProceedToAlgo(elementsList, algo, pbc):
    numAtoms = len(elementsList)
    if algo == 0:
        # How long the ManyMolecule EA takes.
        estTime = round(0.5967 * (numAtoms * numAtoms) - 1.3253 * numAtoms + 3.0362)
    else:
        # How long the PerAtom algorithm takes.
        estTime = round(0.6798 * exp(0.701 * numAtoms))
    sure = messagebox.askquestion(title='Build molecule',
                                  message='Estimated time = {} seconds. Proceed?'.format(estTime))
    if sure == 'yes':
        StartAnalysis(elementsList, algo, pbc)







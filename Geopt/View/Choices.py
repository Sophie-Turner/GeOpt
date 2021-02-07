from tkinter import messagebox
from math import exp
from Controller.Choices import *
from View.Shared import SetUpWindow
from View.Analysis import StartAnalysis

global formula, window, numAtoms


def ChooseFeatures(elementsList, boxText):
    global formula, window, numAtoms
    formula = boxText
    numAtoms = len(elementsList)
    window = tk.Toplevel()
    window.attributes('-topmost', 'true')
    SetUpWindow(window)
    # Show the user their chosen formula.
    lblHelp = Label(window, text=boxText, font=('Agency FB', 18), fg='white', bg='#222222')
    lblHelp.pack(pady=10)
    optionsFrame = tk.Frame(window)

    # The place where options will be chosen.
    Label(optionsFrame, text='Optimisation algorithm', font=('Agency FB', 14)).grid(row=0, column=0)
    # Set up radio buttons.
    # Algorithms for the user to choose from.
    algos = [('Many-molecule evolution', 0, EnableEA),
             ('Per-atom exhaustive test', 1, DisableEA)]
    algo = tk.IntVar()
    # Default MMEA algo.
    algo.set(0)
    for name, value, cmd in algos:
        tk.Radiobutton(optionsFrame, text=name, variable=algo, value=value, command=cmd).grid(row=value+1, column=0)

    # Periodic boundary conditions.
    Label(optionsFrame, text='Periodic boundary conditions', font=('Agency FB', 14)).grid(row=3, column=0)
    pbcs = [('Off', False), ('On', True)]
    pbc = tk.BooleanVar()
    pbc.set(False)
    for name, value in pbcs:
        tk.Radiobutton(optionsFrame, text=name, variable=pbc, value=value).grid(row=value+4, column=0)

    # Info buttons.
    infoBtns = [(lambda: ShowTxtManyMolecule(window), 1, 1), (lambda: ShowTxtPerAtom(window), 2, 1),
                (lambda: ShowTxtPbc(window), 3, 1), (lambda: ShowTxtPopSize(window), 0, 3)]
    for cmd, row, col in infoBtns:
        Button(optionsFrame, text='?', font=('Agency FB bold', 10), command=cmd, bg='yellow').grid(row=row, column=col)

    # Population size chooser.
    defaultEA = tk.BooleanVar()
    popSize = tk.IntVar()
    # The default population size is equal to the number of atoms.
    defaultEA.set(True)
    popSize.set(numAtoms)
    Label(optionsFrame, text='Population size', font=('Agency FB', 14)).grid(row=0, column=2)
    slider = Scale(optionsFrame, from_=2, to=numAtoms + 2, orient=tk.HORIZONTAL, showvalue=0, tickinterval=1,
                   state="disabled", variable=popSize)
    slider.grid(row=2, column=2)
    Checkbutton(optionsFrame, text='Default', variable=defaultEA, onvalue=True, offvalue=False,
                command=lambda: DisableSlider(defaultEA.get(), slider, popSize)).grid(row=1, column=2)

    # Finish buttons.
    Button(optionsFrame, text='Cancel', font=('Agency FB', 14), command=Close).grid(row=5, column=3)
    Button(optionsFrame, text='Build molecule', font=('Agency FB', 14),
           command=(lambda: ProceedToAlgo(elementsList, algo.get(), pbc.get(), popSize.get()))).grid(row=5, column=2)

    optionsFrame.pack()
    window.mainloop()


def DisableSlider(isDefault, slider, popSize):
    # Disable or enable the population size picker.
    if isDefault is True:
        slider['state'] = 'disabled'
        popSize.set(numAtoms)
    else:
        slider['state'] = 'normal'


def Close():
    window.destroy()


def ProceedToAlgo(elementsList, algo, pbc, popSize):
    if algo == 0:
        # How long the ManyMolecule EA takes.
        estTime = round(0.5967 * (numAtoms * numAtoms) - 1.3253 * numAtoms + 3.0362)
    else:
        # How long the PerAtom algorithm takes.
        estTime = round(0.6798 * exp(0.701 * numAtoms))
    sure = messagebox.askquestion(parent=window, title='Build molecule',
                                  message='Estimated time = {} seconds. Proceed?'.format(estTime))
    if sure == 'yes':
        Close()
        StartAnalysis(elementsList, algo, pbc, popSize)







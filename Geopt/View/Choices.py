from tkinter import messagebox
from math import exp
from Controller.Choices import *
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

    # The place where options will be chosen.
    Label(optionsFrame, text='Optimisation algorithm', font=('Agency FB', 14)).grid(row=0, column=0)
    # Set up radio buttons.
    # Algorithms for the user to choose from.
    algos = [('Many-molecule evolution', 0),
             ('Per-atom exhaustive test', 1)]
    algo = tk.IntVar()
    # Default MMEA algo.
    algo.set(0)
    for name, value in algos:
        tk.Radiobutton(optionsFrame, text=name, variable=algo, value=value).grid(row=value+1, column=0)
    # Periodic boundary conditions.
    Label(optionsFrame, text='Periodic boundary conditions', font=('Agency FB', 14)).grid(row=3, column=0)
    pbcs = [('Off', False), ('On', True)]
    pbc = tk.BooleanVar()
    pbc.set(False)
    for name, value in pbcs:
        tk.Radiobutton(optionsFrame, text=name, variable=pbc, value=value).grid(row=value+4, column=0)

    Button(optionsFrame, text='Build molecule', font=('Agency FB', 14),
           command=(lambda: ProceedToAlgo(elementsList, algo.get(), pbc.get()))).grid(row=5, column=2)

    # Info buttons.
    infoBtns = [(lambda: ShowTxtManyMolecule(window), 1), (lambda: ShowTxtPerAtom(window), 2), (lambda: ShowTxtPbc(window), 3)]
    for cmd, row in infoBtns:
        Button(optionsFrame, text='?', font=('Agency FB bold', 10), command=cmd, bg='yellow').grid(row=row, column=1)

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







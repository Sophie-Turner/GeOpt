from math import exp
from Controller.Choices import *
from View.Shared import SetUpWindow
from View.Analysis import StartAnalysis

global formula, window, numAtoms, cores


def ChooseFeatures(elementsList, boxText):
    global formula, window, numAtoms, cores
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

    # Population size chooser.
    defaultEA = tk.BooleanVar()
    popSize = tk.IntVar()
    # The default population size is equal to the number of atoms.
    defaultEA.set(True)
    popSize.set(numAtoms)
    Label(optionsFrame, text='Population size', font=('Agency FB', 14)).grid(row=0, column=2)
    popSlider = Scale(optionsFrame, from_=2, to=numAtoms + 2, orient=tk.HORIZONTAL, showvalue=0, tickinterval=1,
                   state="disabled", variable=popSize)
    popSlider.grid(row=2, column=2)
    Checkbutton(optionsFrame, text='Default', variable=defaultEA, onvalue=True, offvalue=False,
                command=lambda: DisableSlider(defaultEA.get(), popSlider, popSize, 0)).grid(row=1, column=2)

    # Multiprocessing factor chooser.
    findCores = tk.BooleanVar()
    numCores = tk.IntVar()
    findCores.set(True)
    max = 6
    cores = cpu_count()
    if cores is None:
        ShowMessage(window, 6)
        cores = 1
    elif cores > 6:
        max = cores
    numCores.set(cores)
    Label(optionsFrame, text='Parallel processes', font=('Agency FB', 14)).grid(row=3, column=2)
    proSlider = Scale(optionsFrame, from_=1, to=max, orient=tk.HORIZONTAL, showvalue=0, tickinterval=1,
                   state="disabled", variable=numCores)
    proSlider.grid(row=5, column=2)
    Checkbutton(optionsFrame, text='Default', variable=findCores, onvalue=True, offvalue=False,
                command=lambda: DisableSlider(findCores.get(), proSlider, numCores, 1)).grid(row=4, column=2)

    # Graph detail chooser.
    showPosPlot, showPesPlot = tk.BooleanVar(), tk.BooleanVar()
    numPoints = tk.IntVar()
    numPoints.set(300)
    showPlots = [(showPosPlot, 'Show positions tested', 1), (showPesPlot, 'Potential energy surface', 2)]
    Label(optionsFrame, text='Analytical display', font=('Agency FB', 14)).grid(row=0, column=4)
    for plot, text, row in showPlots:
        plot.set(True)
        Checkbutton(optionsFrame, text=text, variable=plot, onvalue=True, offvalue=False,
                    ).grid(row=row, column=4)
    Label(optionsFrame, text='Data points limit').grid(row=3, column=4)
    Scale(optionsFrame, from_=100, to=800, orient=tk.HORIZONTAL, length=200, tickinterval=100,
          variable=numPoints).grid(row=4, column=4, columnspan=2, rowspan=2)

    # Info buttons.
    infoBtns = [(1, 1, 0), (2, 1, 1), (3, 1, 2), (0, 3, 3), (3, 3, 4), (3, 5, 5)]
    for row, col, which in infoBtns:
        Button(optionsFrame, text='?', font=('Agency FB bold', 10), command=lambda which=which: ShowMessage(window, which),
               bg='yellow').grid(row=row, column=col)

    # Finish buttons.
    Button(optionsFrame, text='Build molecule', font=('Agency FB', 14),
           command=(lambda: ProceedToAlgo(elementsList, algo.get(), pbc.get(), popSize.get(), numCores.get(),
                                          showPosPlot.get(), showPesPlot.get(), numPoints.get())))\
        .grid(row=6, column=2)
    Button(optionsFrame, text='Cancel', font=('Agency FB', 14), command=Close).grid(row=6, column=3)

    optionsFrame.pack()
    window.mainloop()


def DisableSlider(isDefault, slider, value, whichSlider):
    # Disable or enable the population size picker.
    values = (numAtoms, cores)
    if isDefault is True:
        slider['state'] = 'disabled'
        # 0 is the population size slider and 1 is the parallel processes slider.
        value.set(values[whichSlider])
    else:
        slider['state'] = 'normal'


def Close():
    window.destroy()


def ProceedToAlgo(elementsList, algo, pbc, popSize, numCores, showPosPlot, showPesPlot, numPoints):
    if algo == 0:
        # How long the ManyMolecule EA takes.
        estTime = round(0.5967 * (numAtoms * numAtoms) - 1.3253 * numAtoms + 3.0362)
    else:
        # How long the PerAtom algorithm takes.
        estTime = round(0.6798 * exp(0.701 * numAtoms))
    sure = messagebox.askquestion(parent=window, title='Build molecule',
                                  message='Estimated time: up to {} seconds. Proceed?'.format(estTime))
    if sure == 'yes':
        Close()
        StartAnalysis(elementsList, algo, pbc, popSize, numCores, showPosPlot, showPesPlot, numPoints)







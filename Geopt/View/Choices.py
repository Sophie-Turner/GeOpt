from math import exp
from ttkwidgets.frames import Balloon
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
    Label(optionsFrame, text='Optimisation algorithm', font=('Agency FB', 14)).grid(row=0, column=0, columnspan=3)
    # Set up radio buttons.
    # Algorithms for the user to choose from.
    algos = [('Many-molecule evolution', 0),
             ('Per-atom exhaustive test', 1)]
    algo = tk.IntVar()
    # Default MMEA algo.
    algo.set(0)
    for name, value in algos:
        tk.Radiobutton(optionsFrame, text=name, variable=algo, value=value, command=lambda: EnableButtons(buttons))\
            .grid(row=value+1, column=0, columnspan=3)

    # Periodic boundary conditions.
    pbcLabel = Label(optionsFrame, text='Periodic boundary conditions', font=('Agency FB', 14))
    pbcLabel.grid(row=3, column=0, columnspan=3, pady=(10, 0))
    #pbcBalloon = Balloon(pbcLabel, headertext=messages[2][0], text=messages[2][1], timeout=0)
    pbcs = [('Off', False), ('On', True)]
    pbc = tk.BooleanVar()
    pbc.set(False)
    for name, value in pbcs:
        tk.Radiobutton(optionsFrame, text=name, variable=pbc, value=value).grid(row=value+4, column=0, columnspan=3)

    # Population size chooser.
    defaultEA = tk.BooleanVar()
    popSize = tk.IntVar()
    # The default population size is equal to the number of atoms.
    defaultEA.set(True)
    popSize.set(numAtoms)
    popLabel = Label(optionsFrame, text='Population size', font=('Agency FB', 14))
    popLabel.grid(row=0, column=4)
    #popBalloon = Balloon(popLabel, headertext=messages[3][0], text=messages[3][1], timeout=0)
    popSlider = Scale(optionsFrame, from_=2, to=numAtoms + 2, orient=tk.HORIZONTAL, showvalue=0, tickinterval=1,
                      variable=popSize)
    popSlider.grid(row=1, column=4)

    # Multiprocessing factor chooser.
    findCores = tk.BooleanVar()
    numCores = tk.IntVar()
    findCores.set(True)
    max = 6
    cores = cpu_count()
    if cores is None:
        ShowMessage(window, -1)
        cores = 1
    elif cores > 6:
        max = cores
    numCores.set(cores)
    proLabel = Label(optionsFrame, text='Parallel processes', font=('Agency FB', 14))
    proLabel.grid(row=2, column=4, pady=(10, 0))
    #proBalloon = Balloon(proLabel, headertext=messages[4][0], text=messages[4][1], timeout=0)
    proSlider = Scale(optionsFrame, from_=1, to=max, orient=tk.HORIZONTAL, showvalue=0, tickinterval=1,
                   state="disabled", variable=numCores)
    proSlider.grid(row=4, column=4)
    Checkbutton(optionsFrame, text='Default', variable=findCores, onvalue=True, offvalue=False,
                command=lambda: DisableSlider(findCores.get(), proSlider, numCores, 1)).grid(row=3, column=4)

    # Graph detail chooser.
    showPosPlot, showPesPlot = tk.BooleanVar(), tk.BooleanVar()
    numPoints = tk.IntVar()
    numPoints.set(300)
    showPlots = [(showPosPlot, 'Show positions tested', 1), (showPesPlot, 'Potential energy surface', 2)]
    Label(optionsFrame, text='Analytical display', font=('Agency FB', 14)).grid(row=0, column=6)
    for plot, text, row in showPlots:
        plot.set(True)
        Checkbutton(optionsFrame, text=text, variable=plot, onvalue=True, offvalue=False,
                    ).grid(row=row, column=6)
    limLabel = Label(optionsFrame, text='Data points limit')
    limLabel.grid(row=3, column=6)
    #limBalloon = Balloon(limLabel, headertext=messages[5][0], text=messages[5][1], timeout=0)
    Scale(optionsFrame, from_=100, to=800, orient=tk.HORIZONTAL, length=200, tickinterval=100,
          variable=numPoints).grid(row=4, column=6, columnspan=2, rowspan=2)

    # Mutation options.
    mutLabel = Label(optionsFrame, text='Radial mutation', font=('Agency FB', 14))
    mutLabel.grid(row=6, column=0, columnspan=3)
    #mutBalloon = Balloon(mutLabel, headertext=messages[6][0], text=messages[6][1], timeout=0)
    # Types of mutation distribution.
    mutations = [('Uniform', 0), ('Gaussian', 1)]
    # Sizes of mutation distributions.
    sizes = [('Very small', -2), ('Small', -1), ('Medium', 0), ('Large', 1), ('Very large', 2)]
    mutDist, mutSize = tk.IntVar(), tk.IntVar()
    # Default mutation distribution type.
    mutDist.set(1)
    mutSize.set(0)
    # Radio buttons.
    for name, value in mutations:
        tk.Radiobutton(optionsFrame, text=name, variable=mutDist, value=value).grid(row=value + 7, column=0, columnspan=3)
    # Mutation size option.
    Label(optionsFrame, text='Distribution size').grid(row=9, column=0, columnspan=3)
    Label(optionsFrame, text='Small').grid(row=10, column=0)
    Label(optionsFrame, text='Large').grid(row=10, column=2)
    sizeSlider = Scale(optionsFrame, from_=-2, to=2, orient=tk.HORIZONTAL, showvalue=0, tickinterval=1, variable=mutSize)
    sizeSlider.grid(row=10, column=1)

    # Permutation option.
    permuteLabel = Label(optionsFrame, text='Permutation', font=('Agency FB', 14))
    permuteLabel.grid(row=5, column=4)
    #permuteBalloon = Balloon(permuteLabel, headertext=messages[7][0], text=messages[7][1], timeout=0)
    permute = tk.BooleanVar()
    permute.set(False)
    permuteOff = tk.Radiobutton(optionsFrame, text='Off', variable=permute, value=False)
    permuteOff.grid(row=6, column=4)
    permuteOn = tk.Radiobutton(optionsFrame, text='On', variable=permute, value=True)
    permuteOn.grid(row=7, column=4)

    # Crossover option.
    crossLabel = Label(optionsFrame, text='Crossover', font=('Agency FB', 14))
    crossLabel.grid(row=8, column=4)
    #crossBalloon = Balloon(crossLabel, headertext=messages[8][0], text=messages[8][1], timeout=0)
    cross = tk.BooleanVar()
    cross.set(False)
    crossOff = tk.Radiobutton(optionsFrame, text='Off', variable=cross, value=False)
    crossOff.grid(row=9, column=4)
    crossOn = tk.Radiobutton(optionsFrame, text='On', variable=cross, value=True)
    crossOn.grid(row=10, column=4)

    # List of buttons whose state depend on algorithm chosen.
    buttons = (permuteOff, permuteOn, crossOff, crossOn)

    # Info buttons.
    infoBtns = [(1, 3, 0), (2, 3, 1)]
    for row, col, which in infoBtns:
        Button(optionsFrame, text='?', font=('Agency FB bold', 10), command=lambda which=which: ShowMessage(window, which),
               bg='yellow').grid(row=row, column=col, padx=(0, 10))

    # Tooltips.
    #tooltips = [pbcLabel, popLabel, proLabel, limLabel, mutLabel, permuteLabel, crossLabel]
    #for i in range(len(tooltips)):
        #Balloon(tooltips[i], headertext=messages[i + 2][0], text=messages[i + 2][1], timeout=0)
    # Finish buttons.
    Button(optionsFrame, text='Restore defaults', font=('Agency FB', 14), command=(lambda: Reset(elementsList, boxText)))\
        .grid(row=9, column=6, rowspan=2)
    Button(optionsFrame, text='Build molecule', font=('Agency FB', 14),
           command=(lambda: ProceedToAlgo(elementsList, algo.get(), pbc.get(), popSize.get(), numCores.get(),
                                          showPosPlot.get(), showPesPlot.get(), numPoints.get(), mutDist.get(),
                                          mutSize.get(), permute.get(), cross.get())))\
        .grid(row=11, column=6, rowspan=2)
    Button(optionsFrame, text='Cancel', font=('Agency FB', 14), command=Close).grid(row=11, column=7, rowspan=2, padx=(0, 10))

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


def EnableButtons(buttons):
    # Disable or enable the EA options.
    firstButton = buttons[0]
    state = 'normal'
    if firstButton['state'] == 'normal':
        state = 'disabled'
    for button in buttons:
        button['state'] = state


def Close():
    window.destroy()


def Reset(elementsList, boxText):
    Close()
    ChooseFeatures(elementsList, boxText)


def ProceedToAlgo(elementsList, algo, pbc, popSize, numCores, showPosPlot, showPesPlot, numPoints, mutDist, mutSize,
                  permute, cross):
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
        try:
            StartAnalysis(elementsList, algo, pbc, popSize, numCores, showPosPlot, showPesPlot, numPoints, mutDist,
                      mutSize, permute, cross)
        except:
            messagebox.showerror(title="Unsupported element", message="Some heavy elements are not yet supported by the energy calculator.\n"
                                                                      "Please try again when the project is finished.")






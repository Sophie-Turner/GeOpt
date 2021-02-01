from View.Shared import SetUpWindow
from Controller.Analysis import *


def ShowInfo(version, pes, refs, atoms, numTests, rank):
    if numTests == 300:
        numTests = '> 300'
    pesData = np.array(pes)
    imageHolders = []

    # Make a larger copy of the PES for a clearer view.
    SurfacePlot(pesData, refs, 3.5, rank)

    window = tk.Toplevel()
    SetUpWindow(window)

    # Title
    Label(window, text='Result #{num} for elements {elements}    Potential energy = {energy} eV'.format(num=rank+1, elements=atoms, energy=version[1]),
          font=('Agency FB', 16), fg='white', bg='#222222').pack(pady=5)

    # Set up grid
    topGrid = Frame(window)
    SetColours(topGrid)
    sizex, sizey = 200, 170

    Label(topGrid, text='Geometric structure', font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
        .grid(row=0, column=0, columnspan=2, padx='5')
    Label(topGrid, text='First 300 positions tested', font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
        .grid(row=0, column=2, columnspan=2, padx='5')

    imageHolders.append(PhotoImage(file="Images/structure{num}.png".format(num=rank)))
    thisImage = imageHolders[-1]
    canvas = tkinter.Canvas(topGrid, width=sizex, height=sizey, bg="#222222")
    canvas.grid(row=1, column=0, padx='5')
    canvas.create_image(sizex / 2, sizey / 2, image=thisImage)

    # Show the stats for this version of the molecule.
    bestInfoBox = Text(topGrid, fg='#EEFFEE', bg="#222222", width="20", height="10")
    bestInfoBox.grid(row=1, column=1, columnspan=1, padx='5')

    imageHolders.append(PhotoImage(file="Images/positions{num}.png".format(num=rank)))
    thisImage = imageHolders[-1]
    canvas = tkinter.Canvas(topGrid, width=sizex, height=sizey)
    canvas.grid(row=1, column=2, padx='5')
    canvas.create_image(sizex / 2, sizey / 2, image=thisImage)

    Label(topGrid, text='Configurations tested: {}'.format(numTests), fg='#EEFFEE', bg='#222222') \
        .grid(row=1, column=3, columnspan=1, padx='5')

    topGrid.pack()

    GetBestInfo(rank, bestInfoBox)

    # Set up bottom grid
    bottomGrid = Frame(window)
    SetColours(bottomGrid)
    sizex, sizey = 400, 340

    Label(bottomGrid, text='Potential energy surfaces found', font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
        .grid(row=0, column=0, columnspan=4, padx='5')

    # Show info for the PES plot.
    pesInfoBox = Text(bottomGrid, fg='#EEFFEE', bg="#222222", width="27", height="18")
    pesInfoBox.grid(row=1, column=0, columnspan=1, padx='5')

    imageHolders.append(PhotoImage(file="Images/pes{num}.png".format(num=rank)))
    thisImage = imageHolders[-1]
    canvas = tkinter.Canvas(bottomGrid, width=sizex, height=sizey)
    canvas.grid(row=1, column=1, padx='5')
    canvas.create_image(sizex/2, sizey/2, image=thisImage)

    # Show a legend for the PES plot.
    legend = Text(bottomGrid, fg='#EEFFEE', bg="#222222", width="6", height="18")
    legend.grid(row=1, column=3, columnspan=1, padx='5')

    bottomGrid.pack()

    SurfaceInfo(pesInfoBox)
    SurfaceLegend(legend, refs)

    window.mainloop()

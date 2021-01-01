from ase.io import write
import tkinter
import matplotlib.pyplot as plt
import numpy as np
from Controller.Analysis import *
from View.Shared import SetUpWindow


def StartAnalysis(elementsList):
    window = tk.Toplevel()
    SetUpWindow(window)
    bestMolecules, population, plot = DoTheEA(elementsList)

    # The largest energy value is needed for scaling.
    eMax = population[0][2]
    firstPlot = plot[0]
    firstPlot = np.array(firstPlot)

    testFig = plt.figure()
    ax = testFig.add_subplot(111, projection='3d')
    # Change the markers to be the symbols and colour according to energy.
    for eachPoint in firstPlot:
        # Calculate energy colours by scaling.
        eCol = eachPoint[0].astype('float64')/eMax
        if eCol > 1:
            eCol = 1
        print('e:', eachPoint[0])
        print('eCol:', eCol)
        ax.scatter(eachPoint[1].astype('float64'), eachPoint[2].astype('float64'), eachPoint[3].astype('float64'),
                   marker='${}$'.format(eachPoint[4]), color=(eCol, 1-eCol, 0.0))
    # Add all 6 versions to the plot.
    plt.savefig("Images/fig3.png")

    for i in range(3):
        fileName = "Images/fig{num}.png".format(num=i)
        write(fileName, bestMolecules[i], rotation='10x,30y,0z')

    picturesFrame = Frame(window)

    photoImages = []
    rows = [2, 4, 4, 2, 1, 3, 3, 1]
    cols = [1, 1, 2, 3]
    sizes = [200, 100, 100, 300]
    colspans = [2, 1, 1, 1]
    rowspans = [1, 1, 1, 3]
    texts = ['Best configuration found   {:.4f} eV'.format(population[0][1]), '#2   {:.4f} eV'.format(population[1][1]),
             '#3   {:.4f} eV'.format(population[2][1]), 'All positions tested']

    for i in range(4):
        size = sizes[i]
        fileName = "Images/fig{num}.png".format(num=i)
        photoImages.append(PhotoImage(file=fileName))
        imageBest = photoImages[i]
        canvas = tkinter.Canvas(picturesFrame, height=size, width=size)
        canvas.grid(row=rows[i], column=cols[i], rowspan=rowspans[i], columnspan=colspans[i], padx='10', pady='5')
        canvas.create_image(size / 2, size / 2, image=imageBest)

    SetColours(picturesFrame)

    for i in range(4):
        Label(picturesFrame, text=texts[i], font=('Agency FB', 16), fg='#EEFFEE', bg='#222222')\
            .grid(row=rows[i+4], column=cols[i], columnspan=colspans[i], pady='5')

    picturesFrame.pack()

    window.mainloop()




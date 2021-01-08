from ase.io import write
import tkinter
import matplotlib.pyplot as plt
import numpy as np
from tkinter import ttk
from Controller.Analysis import *
from View.Shared import SetUpWindow


def PositionPlot(allAtomPlaces, eMax, fileName):
    testFig = plt.figure()
    ax = testFig.add_subplot(111, projection='3d')
    # Change the markers to be the symbols and colour according to energy.
    for eachPoint in allAtomPlaces:
        # Calculate energy colours by scaling. eCol is energy colour.
        eCol = (eachPoint[0].astype('float64') / eMax)
        if eCol > 1:
            eCol = 1
        ax.plot(eachPoint[1].astype('float64'), eachPoint[2].astype('float64'), eachPoint[3].astype('float64'),
                marker='${}$'.format(eachPoint[4]), color=(eCol, 1 - eCol, 0.0))
    plt.savefig("Images/positions{name}.png".format(name=fileName))


def SurfacePlot(pesData, refs, fileName):
    surfFig = plt.figure()
    pax = surfFig.gca(projection='3d')
    for i in range(len(refs)):
        group = refs[i]
        grouped = np.where(pesData[:, 3] == group)
        distances = pesData[grouped, 0].astype('float64')
        energies = pesData[grouped, 1].astype('float64')
        angles = pesData[grouped, 2].astype('float64')
        pax.plot_trisurf(energies[0], distances[0], angles[0])

    pax.view_init(110, -90)
    pax.set_xlabel('Distance between atoms / \u00c5')
    pax.set_ylabel('Potential energy of molecule / eV')
    pax.set_zlabel('Angle of atoms / \u00b0')
    plt.savefig("Images/pes{name}.png".format(name=fileName))


def StartAnalysis(elementsList):
    window = tk.Toplevel()
    SetUpWindow(window)
    bestMolecules, population, plot, pes, refs = DoTheEA(elementsList)

    imageTypes = ['structure', 'pes', 'positions']
    imageHolders = [[], [], []]

    for i in range(3):
        fileName = "Images/structure{num}.png".format(num=i)
        write(fileName, bestMolecules[i], rotation='10x,30y,0z')

        pesData = pes[i]
        pesData = np.array(pesData)
        refData = refs[i]
        SurfacePlot(pesData, refData, i)

        # The largest energy value is needed for scaling.
        eMax = population[i][1]
        allAtomPlaces = plot[i]
        allAtomPlaces = np.array(allAtomPlaces)
        PositionPlot(allAtomPlaces, eMax, i)

    colText = ['Best geometry found\n{:.4f} eV'.format(population[0][2]),
               '2nd best found\n{:.4f} eV'.format(population[1][2]),
               '3rd best found\n{:.4f} eV'.format(population[2][2])]
    rowText = ['', 'Structure', 'Potential energy surfaces found', 'All positions tested']

    gridFrame = Frame(window)
    SetColours(gridFrame)

    # Set up grid
    for i in range(4):
        Label(gridFrame, text=rowText[i], font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
            .grid(row=i, column=0, rowspan=1, columnspan=1, padx='5', pady='5')
        for j in range(3):
            if i == 0:
                Label(gridFrame, text=colText[j], font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
                    .grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5', pady='5')
            else:
                sizex = 200
                sizey = 160
                fileName = "Images/{type}{num}.png".format(num=j, type=imageTypes[i-1])
                imageHolders[i-1].append(PhotoImage(file=fileName))
                thisImage = imageHolders[i-1][j]
                canvas = tkinter.Canvas(gridFrame, height=sizey, width=sizex)
                canvas.grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5', pady='5')
                canvas.create_image(sizex / 2, sizey / 2, image=thisImage)

    gridFrame.pack()

    window.mainloop()

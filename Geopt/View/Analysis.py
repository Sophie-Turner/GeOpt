from ase.io import write
import tkinter
import matplotlib.pyplot as plt
import numpy as np
from Controller.Analysis import *
from View.Shared import SetUpWindow

colours = ['deeppink', 'yellow', 'dodgerblue', 'limegreen', 'darkorange', 'purple', 'red', 'blue']


def PositionPlot(allAtomPlaces, eMax, fileName):
    testFig = plt.figure(figsize=(2.5, 2.5))
    ax = testFig.add_subplot(111, projection='3d')
    # Change the markers to be the symbols and colour according to energy.
    for eachPoint in allAtomPlaces:
        # Calculate energy colours by scaling. eCol is energy colour.
        eCol = (eachPoint[0].astype('float64') / eMax)
        if eCol > 1:
            eCol = 1
        ax.plot(eachPoint[1].astype('float64'), eachPoint[2].astype('float64'), eachPoint[3].astype('float64'),
                marker='${}$'.format(eachPoint[4]), color=(eCol, 1 - eCol, 0.0))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    plt.savefig("Images/positions{name}.png".format(name=fileName))


def SurfacePlot(pesData, refs, eMax, fileName):
    surfFig = plt.figure(figsize=(2.5, 2.5))
    pax = surfFig.gca(projection='3d')
    for i in range(len(refs)):
        if i > 7:
            break
        group = refs[i]
        grouped = np.where(pesData[:, 3] == group)
        distances = pesData[grouped, 0].astype('float64')
        energies = pesData[grouped, 1].astype('float64')
        angles = pesData[grouped, 2].astype('float64')
        pax.plot_trisurf(energies[0], distances[0], angles[0], color=colours[i])

    pax.view_init(110, -90)
    pax.set_xlabel('Distance between atoms / \u00c5', fontsize='x-small')
    pax.set_ylabel('Potential energy of molecule / eV', fontsize='x-small')
    pax.set_zlabel('Angle of atoms / \u00b0', fontsize='x-small')
    axes = [pax.xaxis, pax.yaxis, pax.zaxis]
    for ax in axes:
        for tick in ax.get_major_ticks():
            tick.label.set_fontsize('x-small')
            tick.label.set_rotation('horizontal')

    plt.savefig("Images/pes{name}.png".format(name=fileName))


def StartAnalysis(elementsList):
    window = tk.Toplevel()
    SetUpWindow(window)
    window.geometry("+0+0")
    bestMolecules, population, plot, pes, refs = DoTheEA(elementsList)

    imageTypes = ['structure', 'pes', 'positions']
    imageHolders = [[], [], []]

    for i in range(3):
        fileName = "Images/structure{num}.png".format(num=i)
        write(fileName, bestMolecules[i], rotation='10x,30y,0z')

        # The largest energy value is needed for scaling.
        eMax = population[i][1]

        pesData = pes[i]
        pesData = np.array(pesData)
        refData = refs[i]
        SurfacePlot(pesData, refData, eMax, i)

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
            .grid(row=i, column=0, rowspan=1, columnspan=1, padx='5')
        for j in range(3):
            if i == 0:
                Button(gridFrame, text=colText[j], command=(lambda: SetColours(gridFrame)), font=('Agency FB', 16),
                       fg='#DDFFDD', bg='#555555').grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5')
            else:
                sizex = 300
                sizey = 175
                fileName = "Images/{type}{num}.png".format(num=j, type=imageTypes[i-1])
                imageHolders[i-1].append(PhotoImage(file=fileName))
                thisImage = imageHolders[i-1][j]
                canvas = tkinter.Canvas(gridFrame, height=sizey, width=sizex, bg="#222222")
                canvas.grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5')
                canvas.create_image(sizex / 2, sizey / 2, image=thisImage)
    legend = Text(gridFrame, bg="#222222", width="3", height=len(refData)*1.5)
    legend.grid(row=2, column=4, rowspan=1, columnspan=1, padx='5')
    gridFrame.pack()
    for i in range(len(refData)):
        string = (refData[i])
        legend.tag_configure(string, foreground=colours[i], font=('Agency FB', 14, 'bold'))
        legend.insert(END, string+'\n', string)

    window.mainloop()

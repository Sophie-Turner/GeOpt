from ase.io import write
from Controller.Analysis import *
from View.Shared import SetUpWindow
from View.Info import ShowInfo
from time import time


def StartAnalysis(elementsList, algo, pbc, popSize, numCores, showPosPlot, showPesPlot, numPoints, mutDist, mutSize,
                  permute, cross):
    startTime = time()
    window = tk.Toplevel()
    SetUpWindow(window)
    window.geometry("+0+0")
    bestMolecules, energies, plot, pes, refs = DoTheAlgo(elementsList, algo, pbc, popSize, numCores, numPoints, mutDist,
                                                         mutSize, permute, cross)

    imageTypes = ['structure', 'pes', 'positions']
    if showPesPlot is False:
        imageTypes[1] = 'noPlot'
    if showPosPlot is False:
        imageTypes[2] = 'noPlot'
    imageHolders = [[], [], []]

    if numCores > 2:
        versions = 3
    else:
        versions = numCores
    for i in range(versions):
        fileName = "Images/structure{num}.png".format(num=i)
        write(fileName, bestMolecules[i], rotation='10x,30y,0z')
        fileName = "Images/structure{num}rotated.png".format(num=i)
        write(fileName, bestMolecules[i], rotation='-10x,210y,0z')

        if showPesPlot is True:
            pesData = pes[i]
            pesData = np.array(pesData)
            refData = refs[i]
            SurfacePlot(pesData, refData, 2.6, i)

        allAtomPlaces = ()
        if showPosPlot is True:
            # The largest energy value is needed for scaling.
            eMax = energies[i][0]
            allAtomPlaces = plot[i]
            allAtomPlaces = np.array(allAtomPlaces)
            PositionPlot(allAtomPlaces, eMax, i)

    colText = ['Best geometry found {:.4f} eV\nClick for more info'.format(energies[0][1])]
    if versions > 1:
        colText.append('2nd best found {:.4f} eV\nClick for more info'.format(energies[1][1]))
    if versions > 2:
        colText.append('3rd best found {:.4f} eV\nClick for more info'.format(energies[2][1]))
    rowText = ['', 'Structure', 'Potential energy surfaces found', 'Positions tested']

    gridFrame = Frame(window)
    SetColours(gridFrame)

    # Set up grid
    for i in range(4):
        Label(gridFrame, text=rowText[i], font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
            .grid(row=i, column=0, rowspan=1, columnspan=1, padx='5')
        for j in range(versions):
            if i == 0:
                Button(gridFrame, text=colText[j], command=(lambda j=j: ShowInfo(energies[j], pes[j], refs[j],
                                                                                 elementsList, len(allAtomPlaces), j,
                                                                                 showPosPlot, showPesPlot, numPoints)),
                       font=('Agency FB', 16), fg='#DDFFDD', bg='#555555').grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5')
            else:
                sizex = 300
                sizey = 175
                fileName = "Images/{type}{num}.png".format(num=j, type=imageTypes[i-1])
                imageHolders[i-1].append(PhotoImage(file=fileName))
                thisImage = imageHolders[i-1][j]
                canvas = tk.Canvas(gridFrame, height=sizey, width=sizex, bg="#222222")
                canvas.grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5')
                canvas.create_image(sizex / 2, sizey / 2, image=thisImage)
    # Show a legend for the PES plots.
    if showPesPlot is True:
        legend = Text(gridFrame, fg='#EEFFEE', bg="#222222", width="6", height="11")
        legend.grid(row=2, column=4, rowspan=1, columnspan=1, padx='5')
        gridFrame.pack()
        SurfaceLegend(legend, refData)
    else:
        gridFrame.pack()

    endTime = time()
    print('Time taken = {} seconds'.format(round(endTime - startTime)))

    window.mainloop()

from ase.io import write
from Controller.Analysis import *
from View.Shared import SetUpWindow
from View.Info import ShowInfo
from time import time


def StartAnalysis(elementsList, algo):
    startTime = time()

    window = tk.Toplevel()
    SetUpWindow(window)
    window.geometry("+0+0")
    bestMolecules, energies, plot, pes, refs = DoTheAlgo(elementsList, algo)
    print('length of plot[0] =', len(plot[0]))
    print('length of pes[0] =', len(pes[0]))

    imageTypes = ['structure', 'pes', 'positions']
    imageHolders = [[], [], []]

    for i in range(3):
        fileName = "Images/structure{num}.png".format(num=i)
        write(fileName, bestMolecules[i], rotation='10x,30y,0z')

        # The largest energy value is needed for scaling.
        eMax = energies[i][0]

        pesData = pes[i]
        pesData = np.array(pesData)
        refData = refs[i]
        SurfacePlot(pesData, refData, 2.6, i)

        allAtomPlaces = plot[i]
        allAtomPlaces = np.array(allAtomPlaces)
        PositionPlot(allAtomPlaces, eMax, i)

    colText = ['Best geometry found {:.4f} eV\nClick for more info'.format(energies[0][1]),
               '2nd best found {:.4f} eV\nClick for more info'.format(energies[1][1]),
               '3rd best found {:.4f} eV\nClick for more info'.format(energies[2][1])]
    rowText = ['', 'Structure', 'Potential energy surfaces found', 'Positions tested']

    gridFrame = Frame(window)
    SetColours(gridFrame)

    # Set up grid
    for i in range(4):
        Label(gridFrame, text=rowText[i], font=('Agency FB', 16), fg='#EEFFEE', bg='#222222') \
            .grid(row=i, column=0, rowspan=1, columnspan=1, padx='5')
        for j in range(3):
            if i == 0:
                Button(gridFrame, text=colText[j], command=(lambda j=j: ShowInfo(energies[j], pes[j], refs[j],
                                                                                 elementsList, len(allAtomPlaces), j)),
                       font=('Agency FB', 16), fg='#DDFFDD', bg='#555555').grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5')
            else:
                sizex = 300
                sizey = 175
                fileName = "Images/{type}{num}.png".format(num=j, type=imageTypes[i-1])
                imageHolders[i-1].append(PhotoImage(file=fileName))
                thisImage = imageHolders[i-1][j]
                canvas = tkinter.Canvas(gridFrame, height=sizey, width=sizex, bg="#222222")
                canvas.grid(row=i, column=j+1, rowspan=1, columnspan=1, padx='5')
                canvas.create_image(sizex / 2, sizey / 2, image=thisImage)
    # Show a legend for the PES plots.
    legend = Text(gridFrame, fg='#EEFFEE', bg="#222222", width="6", height="11")
    legend.grid(row=2, column=4, rowspan=1, columnspan=1, padx='5')
    gridFrame.pack()
    SurfaceLegend(legend, refData)

    endTime = time()
    print('Time taken = {} seconds'.format(round(endTime - startTime)))

    window.mainloop()

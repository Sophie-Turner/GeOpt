from ase.io import write
import tkinter
import matplotlib.pyplot as plt
import numpy as np
from Controller.Analysis import *
from View.Shared import SetUpWindow


def StartAnalysis(elementsList):
    window = tk.Toplevel()
    SetUpWindow(window)
    bestMolecules, population, plotE, plotXYZ = DoTheEA(elementsList)

    # This is a placeholder to keep the loop working until an energy plot is made.
    # Create an example graph to practise plotting on canvas...
    x=np.linspace(0,2*np.pi,50)
    y=np.sin(x)
    testFig = plt.figure()
    ax = testFig.add_axes([0,0,1,1])
    ax.plot(x, y)

    for i in range(len(plotE)):
        pass

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
             '#3   {:.4f} eV'.format(population[2][1]), 'Potential energy plot']

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




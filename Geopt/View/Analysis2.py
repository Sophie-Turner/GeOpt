from ase.io import write
from tkinter import *
import tkinter
from PIL import Image, ImageTk
from Controller.Analysis import *
from View.Shared import SetUpWindow


def StartAnalysis(elementsList):
    window = tk.Toplevel()
    SetUpWindow(window)
    molecule, energy = DoTheEA(elementsList)

    write("Images/bestBuildUp.png", molecule, rotation='10x,30y,0z')

    picturesFrame = Frame(window)

    canvases = []
    rows = [2, 4, 4, 2, 1, 3, 3, 1]
    cols = [1, 1, 2, 3]
    sizes = [200, 100, 100, 300]
    colspans = [2, 1, 1, 1]
    rowspans = [1, 1, 1, 3]
    texts = ['Best configuration found   {:.4f} eV'.format(energy),'#2   {:.4f} eV'.format(energy),'#3   {:.4f} eV'.format(energy), 'Potential energy surface']
    imageBest = PhotoImage(file="Images/bestBuildUp.png")

    for i in range(4):
        size = sizes[i]
        canvas = tkinter.Canvas(picturesFrame, height=size, width=size)
        canvas.grid(row=rows[i], column=cols[i], rowspan=rowspans[i], columnspan=colspans[i], padx='10', pady='5')
        canvas.create_image(size/2, size/2, image=imageBest)
        canvases.append(canvas)

    SetColours(picturesFrame)

    for i in range(4):
        Label(picturesFrame, text=texts[i], font=('Agency FB', 16), fg='#EEFFEE', bg='#222222')\
            .grid(row=rows[i+4], column=cols[i], columnspan=colspans[i], pady='5')

    picturesFrame.pack()

    window.mainloop()




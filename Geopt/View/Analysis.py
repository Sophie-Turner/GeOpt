from ase.io import write
from PIL import Image, ImageTk
from tkinter import *
from Controller.Analysis import *
from View.Shared import SetUpWindow

def StartAnalysis(elementsList):
    window = tk.Toplevel()
    SetUpWindow(window)

    molecule, energy = DoTheEA(elementsList)
    pictureFrame = Frame(window)

    write("Images/bestBuildUp.png", molecule, rotation='10x,30y,0z')

    img = PhotoImage(file='Images/bestBuildUp.png')
    Label(pictureFrame, image=img).pack()

    pictureFrame.pack(pady=10)

    window.mainloop()




import tkinter as tk
from tkinter import *
from tkinter import ttk

root = tk.Tk()
root.title('Geopt')
root.resizable(1, 1)

bannerFrame = Frame(root, bg='#222222')
bannerFrame.pack()

lblName = Label(bannerFrame, text='Geopt', font=('Agency FB', 20), fg='limegreen', bg='#222222')
lblName.pack(side=LEFT, padx=20)
lblDescription = Label(bannerFrame, text='Molecular geometry optimisation', font=('Agency FB', 18), fg='white', bg='#222222')
lblDescription.pack(pady=3)

root.mainloop()
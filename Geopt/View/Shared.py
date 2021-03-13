from Controller.Shared import *


def SetUpWindow(window):
    # Get the value of the screen width.
    windowWidth = window.winfo_reqwidth()

    # Stop tkinter putting the window in a ridiculous place.
    positionRight = int(window.winfo_screenwidth() / 4 - windowWidth / 2)
    positionDown = 0
    window.geometry("+{}+{}".format(positionRight, positionDown))

    window.title('Geopt')
    window.resizable(1, 1)
    window.config(bg='#222222')
    bannerFrame = Frame(window, bg='#222222')
    bannerFrame.pack()
    lblName = Label(bannerFrame, text='Geopt', font=('Agency FB', 20), fg='limegreen', bg='#222222')
    lblName.pack(side=LEFT, padx=20)
    lblDescription = Label(bannerFrame, text='Molecular geometry optimisation', font=('Agency FB', 18),
                           fg='white',
                           bg='#222222')
    lblDescription.pack(pady=5)

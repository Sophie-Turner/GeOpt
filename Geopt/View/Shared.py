from Controller.Shared import *

def SetupWindow():
    root = tk.Tk()

    # Get the value of the screen width.
    windowWidth = root.winfo_reqwidth()

    # Stop tkinter putting the window in a ridiculous place.
    positionRight = int(root.winfo_screenwidth() / 4 - windowWidth / 2)
    positionDown = 0
    root.geometry("+{}+{}".format(positionRight, positionDown))

    root.title('Geopt')
    root.resizable(1, 1)
    root.config(bg='#222222')
    bannerFrame = Frame(root, bg='#222222')
    bannerFrame.pack()
    lblName = Label(bannerFrame, text='Geopt', font=('Agency FB', 20), fg='limegreen', bg='#222222')
    lblName.pack(side=LEFT, padx=20)
    lblDescription = Label(bannerFrame, text='Molecular geometry optimisation', font=('Agency FB', 18), fg='white', bg='#222222')
    lblDescription.pack(pady=3)
    return root
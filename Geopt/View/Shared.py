from Controller.Shared import *


class Application:
    # Don't let the application be duplicated.
    __instance = None

    def __init__(self):
        if Application.__instance is None:
            Application.__instance = self
            self.root = tk.Tk()

            # Get the value of the screen width.
            windowWidth = self.root.winfo_reqwidth()

            # Stop tkinter putting the window in a ridiculous place.
            positionRight = int(self.root.winfo_screenwidth() / 4 - windowWidth / 2)
            positionDown = 0
            self.root.geometry("+{}+{}".format(positionRight, positionDown))

            self.root.title('Geopt')
            self.root.resizable(1, 1)
            self.root.config(bg='#222222')
            bannerFrame = Frame(self.root, bg='#222222')
            bannerFrame.pack()
            lblName = Label(bannerFrame, text='Geopt', font=('Agency FB', 20), fg='limegreen', bg='#222222')
            lblName.pack(side=LEFT, padx=20)
            lblDescription = Label(bannerFrame, text='Molecular geometry optimisation', font=('Agency FB', 18), fg='white',
                                   bg='#222222')
            lblDescription.pack(pady=3)

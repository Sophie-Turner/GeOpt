from Controller.Shared import *

def SetupWindow():
    root = tk.Tk()
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
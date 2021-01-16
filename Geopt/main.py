import tkinter as tk
from tkinter import *
from View.BuildNew import MakeUI
from View.Shared import SetUpWindow


class Application:
    # Don't let the application be duplicated.
    __instance = None

    def __init__(self):
        if Application.__instance is None:
            Application.__instance = self
            self.root = tk.Tk()


if __name__ == '__main__':
    application = Application()
    root = application.root
    MakeUI(root)

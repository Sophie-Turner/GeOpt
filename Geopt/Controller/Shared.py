import tkinter as tk
from tkinter import *

def SetColours(frame):
    frame.config(bg='#222222')
    for item in frame.children.values():
        item.config(bg='#555555', fg='#DDFFDD')
# Functions associated with the BuildNew View
import xml.etree.ElementTree as ET
import tkinter as tk
from tkinter import *
from tkinter import simpledialog, messagebox
from Model.Molecules import Molecule
from Controller.Shared import *

elementList = []
formula = []
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

def GetXML():
    treeMain = ET.parse('../Model/mainblocks.xml')
    treeF = ET.parse('../Model/fblock.xml')
    treerootMain = treeMain.getroot()
    treerootF = treeF.getroot()
    return treerootMain, treerootF


def AddElement(frame, element):
    # Build up our molecule from elements
    howMany = 1
    name = element[1].text
    if(name=='Hydrogen' or name=='Carbon' or name=='Oxygen'):
        howMany = simpledialog.askinteger(title=name, prompt='How many {} atoms?'.format(name))
        if not howMany:
            return

    if howMany >= 20:
        messagebox.showerror(title="Error", message="Too many {}s!".format(name))
        return
    for eachatom in range(howMany):
        elementList.append(element)
    # Update the displayed info. Separated into another function for readability
    UpdateLabels(frame, element, howMany, name)
    frame.config(bg='#222222')
    frame.pack()


def UpdateLabels(frame, element, howMany, name):
    # Remove existing labels
    ClearLabels(frame)

    # Put new formula label in
    symbol = element[0].text
    if howMany == 1:
        formula.append(symbol)
    else:
        thisString = symbol + str(howMany).translate(subscript)
        formula.append(thisString)
    lblMolecule = Label(frame, text=formula, font=('Agency FB', 30))
    lblMolecule.config(bg='#222222', fg='#EEFFEE')
    lblMolecule.pack()


def Build(box, isEntry):
    # If they chose elements from the table
    if isEntry is False:
        if elementList:
            thisMolecule = Molecule(elementList, None, None, None)
        else:
            messagebox.showinfo(title='No elements',
                                message='Please choose the elements for the molecule.')
    # If they used the entry box instead
    elif isEntry is True:
        if len(box.get()) == 0:
            # If the user hasn't entered any text, see if they selected elements from the table
            isEntry = False
            Build(box, isEntry)
        else:
            # build up element list from entry box
            pass


def Clear(frame, box):
    elementList = []
    formula = []
    box.delete(0, END)
    ClearLabels(frame)


def ClearLabels(frame):
    try:
        for label in frame.children.values():
            label.destroy()
            ClearLabels(frame)
    except:
        frame.pack_forget()


# how to get certain XML elements:
    #for element in treeroot:
        #print(element.tag, element.attrib)
        #print(element[2].text)
    #print(root[0][2].text)
    #print("len(treeroot):", len(treeroot))
    #print("len(treeroot[0]):", len(treeroot[0]))
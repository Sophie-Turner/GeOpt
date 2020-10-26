# Functions associated with the BuildNew View
import xml.etree.ElementTree as ET
import tkinter as tk
from tkinter import *
from tkinter import simpledialog, messagebox
from Model.Molecules import Molecule
from Controller.Shared import *

elementList = []

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
    for eachatom in range(howMany):
        elementList.append(element)
    if howMany == 1:
        lblAdded = Label(frame, text="{} was added to the molecule.".format(name),
                        font=('Agency FB', 14))
        lblAdded.config(bg='#222222', fg='#EEFFEE')
        lblAdded.pack()
        frame.pack()
    elif howMany < 20:
        lblAdded = Label(frame, text="{} {}s were added to the molecule.".format(howMany, name),
                        font=('Agency FB', 14))
        lblAdded.config(bg='#222222', fg='#EEFFEE')
        lblAdded.pack()
        frame.pack()
    else:
        messagebox.showerror(title="Error", message="Too many {}s!".format(name))
    frame.config(bg='#222222')



def Build():
    if elementList:
        thisMolecule = Molecule(elementList, None, None, None)


def Clear(frame, box):
    elementList = []
    box.delete(0, END)
    try:
        for label in frame.children.values():
            label.destroy()
            Clear(frame, box)
    except:
        frame.pack_forget()


# how to get certain XML elements:
    #for element in treeroot:
        #print(element.tag, element.attrib)
        #print(element[2].text)
    #print(root[0][2].text)
    #print("len(treeroot):", len(treeroot))
    #print("len(treeroot[0]):", len(treeroot[0]))
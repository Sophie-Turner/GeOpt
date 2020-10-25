# Functions associated with the BuildNew View
import xml.etree.ElementTree as ET
import tkinter as tk
from tkinter import *
from tkinter import messagebox

elementList = []

def GetXML():
    treeMain = ET.parse('../Model/mainblocks.xml')
    treeF = ET.parse('../Model/fblock.xml')
    treerootMain = treeMain.getroot()
    treerootF = treeF.getroot()
    return treerootMain, treerootF


def AddElement(frame, element):
    # Build up our molecule from elements
    if(element[0]=='H' or element[0]=='C'):
        pass
    try:
        elementList.append(element)
        lblName = Label(frame, text="{} was added to the molecule.".format(element[1].text),
                        font=('Agency FB', 14))
        lblName.pack()
    except:
        messagebox.showerror(title="Error",
                             message="{} could not be added to the molecule.".format(element[1].text))



# how to get certain XML elements:
    #for element in treeroot:
        #print(element.tag, element.attrib)
        #print(element[2].text)
    #print(root[0][2].text)
    #print("len(treeroot):", len(treeroot))
    #print("len(treeroot[0]):", len(treeroot[0]))
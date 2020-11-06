# Functions associated with the BuildNew View
import xml.etree.ElementTree as ET
from tkinter import simpledialog, messagebox
from Model.Molecules import Molecule
from Controller.Shared import *

formula = []
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
backToNormal = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def GetXML():
    treeMain = ET.parse('../Model/mainblocks.xml')
    treeF = ET.parse('../Model/fblock.xml')
    treerootMain = treeMain.getroot()
    treerootF = treeF.getroot()
    return treerootMain, treerootF


def AddElement(box, element):
    # Build up our molecule from elements
    howMany = 1
    name = element[1].text
    if name=='Hydrogen' or name=='Carbon' or name=='Oxygen':
        howMany = simpledialog.askinteger(title=name, prompt='How many {} atoms?'.format(name))
        if not howMany:
            return
    if howMany > 20:
        messagebox.showerror(title="Error", message="Too many {}s!".format(name))
        return
    # Update the displayed info. Separated into another function for readability
    UpdateFormula(box, element, howMany, name)


def UpdateFormula(box, element, howMany, name):

    # Put new formula in
    symbol = element[0].text
    if howMany == 1:
        formula.append(symbol)
    else:
        thisString = symbol + str(howMany).translate(subscript)
        formula.append(thisString)

    box.delete(0, END)
    box.insert(END, formula)


def Build(box):
    boxText = box.get()
    boxLength = len(boxText)
    # Remove all white spaces
    boxText = boxText.replace(" ", "")
    if boxLength == 0:
        messagebox.showinfo(title='No elements',
                            message='Please choose the elements for the molecule.')
    else:
        # Disassemble the string and piece it back together one atom at a time
        elementsList = []
        thisAtom = ''
        for character in boxText:
            if character.isupper():
                elementsList.append(thisAtom)
                thisAtom = character
            elif character.isdigit():
                character = character.translate(backToNormal)
                howMany = int(character)
                for i in range(howMany-1):
                    elementsList.append(thisAtom)
            else:
                thisAtom = thisAtom + character
        elementsList.append(thisAtom)
        del elementsList[0]
        print(elementsList)

        # Create an instance of the molecule class and create the molecule.
        thisMolecule = Molecule(elementsList, None, None, None)
        thisMolecule.ModelMolecule()




def Clear(box):
    formula.clear()
    box.delete(0, END)


# How to remove all labels from a frame
# def ClearLabels(frame):
#     try:
#         for label in frame.children.values():
#             label.destroy()
#             ClearLabels(frame)
#     except:
#         frame.pack_forget()


# How to get certain XML elements
    #for element in treeroot:
        #print(element.tag, element.attrib)
        #print(element[2].text)
    #print(root[0][2].text)
    #print("len(treeroot):", len(treeroot))
    #print("len(treeroot[0]):", len(treeroot[0]))
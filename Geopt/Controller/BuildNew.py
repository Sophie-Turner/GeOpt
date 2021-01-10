# Functions associated with the BuildNew View
from tkinter import simpledialog, messagebox
from Model.InteractWithData import GetXML
from Controller.Shared import *
from View.Analysis import StartAnalysis

formula = []
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
backToNormal = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def AddElement(box, element):
    # Build up our molecule from elements
    howMany = 1
    name = element[1].text
    if name=='Hydrogen' or name=='Carbon' or name=='Oxygen':
        howMany = simpledialog.askinteger(title=name, prompt='How many {} atoms?'.format(name))
        if not howMany:
            return
    if howMany > 49:
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
        howMany = 0
        for i in range(len(boxText)):
            # Don't let the user make huge molecules!
            if len(elementsList) > 20:
                messagebox.showerror(title='Large molecule',
                                     message='This molecule is too large to build! Cancelling...')
                Clear(box)
                return
            # Skip a loop iteration if the last character was 2-digit to prevent overwriting.
            if howMany > 9:
                howMany = 0
            else:
                # We'll need to check one node ahead.
                character = boxText[i]
                try:
                    nextCharacter = boxText[i+1]
                except:
                    nextCharacter = ''
                # Capital letter is the start of a new element so add the previous one to the list.
                if character.isupper():
                    elementsList.append(thisAtom)
                    thisAtom = character
                # Check for multiples.
                elif character.isdigit():
                    character = character.translate(backToNormal)
                    nextCharacter = nextCharacter.translate(backToNormal)
                    howMany = int(character)
                    # It might be a two-digit number.
                    if nextCharacter.isdigit():
                        howMany = howMany * 10 + int(nextCharacter)
                    for j in range(howMany-1):
                        elementsList.append(thisAtom)
                else:
                    # It could be an element with two letters in the symbol.
                    thisAtom = thisAtom + character
        elementsList.append(thisAtom)
        del elementsList[0]

        #try:
        StartAnalysis(elementsList)
        #except:
            #messagebox.showerror(title="Invalid input", message="Please enter a valid molecular formula, e.g. H2SO4")
            #elementsList.clear()
            #Clear(box)


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



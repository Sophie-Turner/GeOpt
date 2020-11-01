from View.Shared import *
from Controller.BuildNew import *
from Controller.Shared import *

root = SetupWindow()

def MakeTable(tableFrame, periods, xmlList, textFrame):
    for i in range(periods):
        for j in xmlList[i]:
            symbol = j[0].text
            group = j[4].text
            element = Button(tableFrame, text=symbol, command=lambda j=j: AddElement(textFrame, j), font=('Agency FB', 12))
            element.grid(row=i, column=group)
            element.config(width=5)
    SetColours(tableFrame)


# get list from XML
elements = GetXML()
mainBlock = elements[0]
fBlock = elements[1]

lblHelp = Label(root, text='Select elements or type the molecular formula.',
                       font=('Agency FB', 16), fg='white', bg='#222222')
lblHelp.pack(pady=10)

formulaFrame = Frame(root)
lblFormula = Label(formulaFrame, text='Molecular formula:', font=('Agency FB', 14))
lblFormula.grid(row=1, column=1, padx='5')
entryFormula = Entry(formulaFrame, font=('Agency FB', 14))
entryFormula.grid(row=1, column=2, padx='5')
btnBuild = Button(formulaFrame, text='Build molecule', command=lambda: Build(entryFormula, True), font=('Agency FB', 14))
btnBuild.grid(row=1, column=3, padx='5')
SetColours(formulaFrame)
lblFormula.config(bg='#222222', fg='#EEFFEE')
formulaFrame.pack()

tableFrame1 = Frame(root)
tableFrame1.pack(pady=10)

tableFrame2 = Frame(root)
tableFrame2.pack(pady=10)

buildFrame = Frame(root)
btnBuild = Button(buildFrame, text='Build molecule', command=lambda: Build(entryFormula, False), font=('Agency FB', 16))
btnBuild.grid(row=1, column=3, padx='5')
btnClear = Button(buildFrame, text='Clear', command=lambda: Clear(textFrame, entryFormula), font=('Agency FB', 16))
btnClear.grid(row=1, column=4, padx='5')
SetColours(buildFrame)
buildFrame.pack()

textFrame = Frame(root)
textFrame.pack()

# populate periodic table
periods = 7
MakeTable(tableFrame1, periods, mainBlock, textFrame)
periods = 2
MakeTable(tableFrame2, periods, fBlock, textFrame)


root.mainloop()



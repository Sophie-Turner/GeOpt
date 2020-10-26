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

tableFrame1 = Frame(root)
tableFrame1.pack(pady=10)

tableFrame2 = Frame(root)
tableFrame2.pack(pady=10)

buildFrame = Frame(root)
entryFormula = Entry(buildFrame, font=('Agency FB', 14))
entryFormula.grid(row=1, column=2, padx='5')
btnBuild = Button(buildFrame, text='Build molecule', command=Build, font=('Agency FB', 16))
btnBuild.grid(row=1, column=3, padx='5')
btnClear = Button(buildFrame, text='Clear', command=lambda: Clear(textFrame, entryFormula), font=('Agency FB', 16))
btnClear.grid(row=1, column=4, padx='5')
SetColours(buildFrame)
lblFormula = Label(buildFrame, text='Molecular formula:', font=('Agency FB', 14))
lblFormula.config(bg='#222222', fg='#EEFFEE')
lblFormula.grid(row=1, column=1, padx='5')
buildFrame.pack()

textFrame = Frame(root)
textFrame.pack()

# populate periodic table
periods = 7
MakeTable(tableFrame1, periods, mainBlock, textFrame)
periods = 2
MakeTable(tableFrame2, periods, fBlock, textFrame)


root.mainloop()



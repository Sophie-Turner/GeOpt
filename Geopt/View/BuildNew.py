from View.Shared import *
from Controller.BuildNew import *

root = SetupWindow()

def MakeTable(frame, periods, xmlList):
    for i in range(periods):
        for j in xmlList[i]:
            symbol = j[0].text
            group = j[4].text
            element = Button(frame, text=symbol, command=lambda j=j: AddElement(root, j), font=('Agency FB', 12))
            element.grid(row=i, column=group)
            element.config(width=5)



# get list from XML
elements = GetXML()
mainBlock = elements[0]
fBlock = elements[1]

# Populate the main blocks
periods = 7
tableFrame1 = Frame(root)
tableFrame1.pack(pady=10)
MakeTable(tableFrame1, periods, mainBlock)

# Populate the f block
periods = 2
tableFrame2 = Frame(root)
tableFrame2.pack(pady=10)
MakeTable(tableFrame2, periods, fBlock)

root.mainloop()



from View.Shared import *
from Controller.BuildNew import *


def MakeTable(frame, periods, xmlList):
    for i in range(periods):
        for j in xmlList[i]:

            symbol = j[0].text
            group = j[4].text
            element = Button(frame, text=symbol, command=TestClick, font=('Agency FB', 12))
            element.grid(row=i, column=group)
            element.config(width=5)


def TestClick():
    print("Clicked the button.")


# get list from XML
list = GetXML()

root = SetupWindow()

tableFrame = Frame(root)
tableFrame.pack()

# Populate the main blocks
periods = 4
MakeTable(tableFrame, periods, list)

# Populate the f block
# We will need to change the way it finds them as this will find the first 2 periods atm, not the last 2
#periods = 2
#MakeTable(tableFrame, periods, list)

root.mainloop()



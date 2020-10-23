from View.Shared import *
from Controller.BuildNew import *


def MakeTable(frame, periods, groups, xmlList):
    for i in range(periods):
        for j in range(groups):
            element = Button(frame, text=xmlList[i][j], command=TestClick, font=('Agency FB', 12))
            element.grid(row=i, column=j)
            element.config(width=5)


def TestClick():
    print("Clicked the button.")


# get list from XML
# just testing it with this for now...
list = [("A", "", "", "B"),
        ("C", "D", "E", "F"),
        ("G", "H", "I", "J")]

groups = len(list[0])
periods = len(list)

root = SetupWindow()

tableFrame = Frame(root)
tableFrame.pack()

MakeTable(tableFrame, periods, groups, list)

root.mainloop()



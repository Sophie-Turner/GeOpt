# functions associated with the analysis view
from Controller.Shared import *
# from Model.EAmanyMolecules import StartEA
from Model.EAperAtom import StartEA


def DoTheEA(elementsList):
    bestMolecules, population, plot = StartEA(elementsList)
    return bestMolecules, population, plot

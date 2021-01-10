# functions associated with the analysis view
from Controller.Shared import *
# from Model.EAmanyMolecules import StartEA
from Model.EAperAtom import StartEA
import tkinter
import matplotlib.pyplot as plt
import numpy as np

colours = ['deeppink', 'yellow', 'dodgerblue', 'limegreen', 'darkorange', 'purple', 'red', 'blue']


def DoTheEA(elementsList):
    bestMolecules, population, plot, pes, refs = StartEA(elementsList)
    return bestMolecules, population, plot, pes, refs


def PositionPlot(allAtomPlaces, eMax, fileName):
    testFig = plt.figure(figsize=(2.7, 2.7))
    ax = testFig.add_subplot(111, projection='3d')
    # Change the markers to be the symbols and colour according to energy.
    for eachPoint in allAtomPlaces:
        # Calculate energy colours by scaling. eCol is energy colour.
        eCol = (eachPoint[0].astype('float64') / eMax)
        if eCol > 1:
            eCol = 1
        ax.plot(eachPoint[1].astype('float64'), eachPoint[2].astype('float64'), eachPoint[3].astype('float64'),
                marker='${}$'.format(eachPoint[4]), color=(eCol, 1 - eCol, 0.0))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    plt.savefig("Images/positions{name}.png".format(name=fileName))


def SurfacePlot(pesData, refs, size, fileName):
    surfFig = plt.figure(figsize=(size, size))
    pax = surfFig.gca(projection='3d')
    for i in range(len(refs)):
        if i > 7:
            break
        group = refs[i]
        grouped = np.where(pesData[:, 3] == group)
        distances = pesData[grouped, 0].astype('float64')
        energies = pesData[grouped, 1].astype('float64')
        angles = pesData[grouped, 2].astype('float64')
        pax.plot_trisurf(energies[0], distances[0], angles[0], color=colours[i])

    pax.view_init(110, -90)
    pax.set_xlabel('Distance between atoms / \u00c5', fontsize='x-small')
    pax.set_ylabel('Potential energy of molecule / eV', fontsize='x-small')
    pax.set_zlabel('Angle of atoms / \u00b0', fontsize='x-small')
    axes = [pax.xaxis, pax.yaxis, pax.zaxis]
    for ax in axes:
        for tick in ax.get_major_ticks():
            tick.label.set_fontsize('x-small')
            tick.label.set_rotation('horizontal')

    plt.savefig("Images/pes{name}.png".format(name=fileName))

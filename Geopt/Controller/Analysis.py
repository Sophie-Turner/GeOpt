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
    global surfData, surfRefs
    surfData, surfRefs = pesData, refs
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
        pax.plot_trisurf(distances[0], energies[0], angles[0], color=colours[i])

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


def SurfaceInfo(infoBox):
    infoBox.insert(END, "Energy minima:\n")
    for i in range(len(surfRefs)):
        if i > 7:
            break
        group = surfRefs[i]
        grouped = np.where(surfData[:, 3] == group)
        distances = (surfData[grouped, 0].astype('float64'))[0]
        energies = (surfData[grouped, 1].astype('float64'))[0]
        angles = (surfData[grouped, 2].astype('float64'))[0]
        minE = min(energies)
        minI = np.where(energies == minE)
        atomsDistance = group[0:len(group)]
        distance = str(distances[minI])
        angle = str(angles[minI])
        txt = "Distance between {} {}\u00c5\nAngle over {} {}\u00b0\nPotential energy {} eV\n\n".format\
            (atomsDistance, distance, group, angle, minE)
        txt = txt.replace("]", "")
        txt = txt.replace("[", "")
        infoBox.insert(END, txt)


def SurfaceLegend(legendBox, refData):
    for i in range(len(refData)):
        string = (refData[i])
        legendBox.tag_configure(string, foreground=colours[i], font=('Agency FB', 14, 'bold'))
        legendBox.insert(END, string + '\n', string)
        if i == 7:
            break

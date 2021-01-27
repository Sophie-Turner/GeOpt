# functions associated with the analysis view
from Controller.Shared import *
from Model.EAmanyMolecules import StartEA
from Model.PerAtom import Start
import tkinter
import matplotlib.pyplot as plt
import numpy as np

colours = ['deeppink', 'yellow', 'dodgerblue', 'limegreen', 'darkorange', 'purple', 'red', 'blue']
global bestVersions, surfData, surfRefs


def DoTheAlgo(elementsList):
    plt.close('all')
    global bestVersions
    bestVersions, energies, plot, pes, refs = StartEA(elementsList)
    return bestVersions, energies, plot, pes, refs


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
        if i < 8:
            distances, energies, angles, group = GetGroup(i)
            pax.plot_trisurf(distances, energies, angles, color=colours[i])
        else:
            break

    pax.view_init(110, -90)
    pax.set_xlabel('Distance between atoms / \u00c5', fontsize='x-small')
    pax.set_ylabel('Potential energy of molecule / eV', fontsize='x-small')
    pax.set_zlabel('Angle / \u00b0', fontsize='x-small')
    axes = [pax.xaxis, pax.yaxis, pax.zaxis]
    for ax in axes:
        for tick in ax.get_major_ticks():
            tick.label.set_fontsize('x-small')
            tick.label.set_rotation('horizontal')

    plt.savefig("Images/pes{name}.png".format(name=fileName))


def SurfaceInfo(infoBox):
    infoBox.insert(END, "Interaction energy minima:\n\n")
    for i in range(len(surfRefs)):
        if i < 8:
            distances, energies, angles, group = GetGroup(i)
            minE = min(energies)
            minI = np.where(energies == minE)
            atom1, atom2 = group[0], group[1]
            distance = str(distances[minI])
            angle = str(angles[minI])
            txt = "Distance between {}{} {} \u00c5\nAngle over {} {}\u00b0\nPotential energy {} eV\n\n".format\
                (atom1, atom2, distance, group, angle, minE)
            txt = txt.replace("]", "")
            txt = txt.replace("[", "")
            infoBox.insert(END, txt)
        else:
            break
    infoBox.configure(state="disabled")


def SurfaceLegend(legendBox, refData):
    legendBox.insert(END, 'Legend')
    for i in range(len(refData)):
        if i < 8:
            string = (refData[i])
            legendBox.tag_configure(string, foreground=colours[i], font=('Agency FB', 14, 'bold'))
            legendBox.insert(END, string + '\n', string)
        else:
            break
    legendBox.configure(state="disabled")


def GetBestInfo(rank, txtBox):
    # Find distances and angles.
    # Populate box in the same loops as testing properties to make it faster.
    molecule = bestVersions[rank]
    refAngs, refDists = [], []
    numAtoms = len(molecule)
    txtBox.insert(END, "This shape:\n\n")
    for i in range(numAtoms):
        for j in range(numAtoms):
            if i != j:
                atom1, atom2 = molecule[i].symbol, molecule[j].symbol
                refDist = atom1 + atom2
                # Don't repeat the same information.
                # If there is a distance, show it.
                if refDist not in refDists:
                    refDists.append(refDist)
                    distance = molecule.get_distance(i, j)
                    txt = "Distance between {}{} {} \u00c5\n\n".format(atom1, atom2, distance)
                    txt = txt.replace("]", "")
                    txt = txt.replace("[", "")
                    txtBox.insert(END, txt)

                # If there's an angle, show it.
                if j < numAtoms-1:
                    jNext = j+1
                    if jNext != i:
                        atom3 = molecule[jNext].symbol
                        refAng = atom1 + atom2 + atom3
                        # Don't repeat the same information.
                        if refAng not in refAngs:
                            refAngs.append(refAng)
                            angle = molecule.get_angle(i, j, jNext)
                            txt = "Angle over {}{}{} {}\u00b0\n\n".format(atom1, atom2, atom3, angle)
                            txt = txt.replace("]", "")
                            txt = txt.replace("[", "")
                            txtBox.insert(END, txt)
    txtBox.configure(state="disabled")


def GetGroup(n):
    group = surfRefs[n]
    grouped = np.where(surfData[:, 3] == group)
    distances = (surfData[grouped, 0].astype('float64'))[0]
    energies = (surfData[grouped, 1].astype('float64'))[0]
    angles = (surfData[grouped, 2].astype('float64'))[0]
    return distances, energies, angles, group


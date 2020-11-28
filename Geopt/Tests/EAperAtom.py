from concurrent import futures
from Tests.EAs import *

global arrayChanges
global energyCalculations

def StartEA(elementsList):
    # Set up initial values & placeholders
    # calc = SetUpVasp()
    calc = EMT()
    overallBestEnergy = 1000
    bestMolecule = None

    # Get the final box size
    thisPopulation = Population(elementsList)
    boxSize = thisPopulation.boxSize
    length = boxSize[0]
    gaussRange = [length/16, length/8, length/4, length/len(elementsList)]

    if __name__ == '__main__':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, gaussRange[i], calc) for i in range(4)]
            for f in futures.as_completed(results):
                thisResult = f.result()
                thisEnergy = thisResult[1]
                if thisEnergy < overallBestEnergy:
                    overallBestEnergy = thisEnergy
                    bestMolecule = thisResult[0]

        newMolecule = Atoms(bestMolecule, cell=boxSize)
        newMolecule.center()
        write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/bestBuildUp.png", newMolecule,
            rotation='10x,30y,0z')
        print("newMolecule", newMolecule)
        print("best energy:", overallBestEnergy)


def Evolve(elementsList, boxSize, gaussRange, calc):
    buildUp = []
    firstCoordinates = [0, 0, 0]
    for eachAtom in elementsList:
        # Add each atom one at a time to the molecule.
        atomsObject = Atom(eachAtom, firstCoordinates)
        # Keep track of the additions.
        buildUp.append(atomsObject)
        # Reset the energy test as it is likely to increase as more atoms are added.
        buildUpBestEnergy = 1000
        # Move every atom around every other atom.
        for eachAtomToMove in buildUp:
            for eachAtomFixed in buildUp:
                for i in range(4):
                    MoveOneAtom(eachAtomFixed, eachAtomToMove, gaussRange)
                    newMolecule = Atoms(buildUp, cell=boxSize)
                    newMolecule.center()
                    currentEnergy = GetEnergy(newMolecule, calc)
                    if currentEnergy < buildUpBestEnergy:
                        buildUpBestEnergy = currentEnergy
                        buildUpBestStructure = buildUp
        buildUp = buildUpBestStructure
        newMolecule = Atoms(buildUp, cell=boxSize)
        newMolecule.center()
    return buildUp, buildUpBestEnergy


def MoveOneAtom(fixedAtom, atomToMove, sigma):
    mean = fixedAtom.position
    atomToMove.position = (random.normal(mean, sigma, 3))


StartEA(testWater)
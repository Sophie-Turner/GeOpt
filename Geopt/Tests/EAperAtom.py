from Tests.EAs import *

def StartEA(elementsList):
    # Set up initial values & placeholders
    # calc = SetUpVasp()
    calc = EMT()
    overallBestEnergy = 1000
    bestMolecule = None

    # Move Hydrogens to the end of the list, as they can mess up the calculations at the start.
    moves = 0
    while moves < len(elementsList):
        for eachAtom in elementsList:
            if eachAtom == 'H':
                elementsList.remove(eachAtom)
                elementsList.append(eachAtom)
        moves += 1

    # Get the final box size
    thisPopulation = Population(elementsList)
    boxSize = thisPopulation.boxSize
    covRads = thisPopulation.covRads

    if __name__ == '__main__':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, covRads, calc) for i in range(1)]
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
        print("best energy:", overallBestEnergy)


def Evolve(elementsList, boxSize, covRads, calc):
    buildUp = []
    buildUpBestStructure = buildUp
    firstCoordinates = [0, 0, 0]
    onlyH = True
    for eachAtom in elementsList:
        # Add each atom one at a time to the molecule.
        atomsObject = Atom(eachAtom, firstCoordinates)
        # Keep track of the additions.
        buildUp.append(atomsObject)
        # See if buildUp ONLY contains H, as this is a special case.
        for atom in buildUp:
            if atom.symbol != 'H':
                onlyH = False
        # Reset the energy test as it is likely to increase as more atoms are added.
        buildUpBestEnergy = 1000
        # Move each atom around every atom.
        buildUpBestEnergy, buildUpBestStructure = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure,
                                                                onlyH, covRads, boxSize, calc)
        buildUp = buildUpBestStructure
    print("this best energy:", buildUpBestEnergy)

    similarity = 0
    iterations = 1
    while similarity < 2:
        # Move each atom around every atom.
        newBestEnergy, newBestStructure = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure,
                                                                onlyH, covRads, boxSize, calc)
        if abs(newBestEnergy - buildUpBestEnergy) < 0.5:
            similarity += 1
        else:
            similarity = 0
        if newBestEnergy < buildUpBestEnergy:
            buildUpBestEnergy = newBestEnergy
            buildUpBestStructure = newBestStructure
            buildUp = newBestStructure
        iterations += 1
        print("this best energy:", newBestEnergy)
    print("iterations:", iterations)
    return buildUp, buildUpBestEnergy


def TestAllPlaces(buildUp, bestEnergy, bestStructure, onlyH, covRads, boxSize, calc):
    # Move each atom around every atom.
    for eachAtomToMove in range(len(buildUp)):
        for eachAtomFixed in range(len(buildUp)):
            fixed = buildUp[eachAtomFixed]
            move = buildUp[eachAtomToMove]
            # Don't let H2 form and break the molecule apart
            # but if the molecule only contains H (eg. H2) we still need to do it.
            if fixed.symbol != 'H' or onlyH is True:
                moveRange = covRads[eachAtomToMove] + covRads[eachAtomFixed]
                # Multithread this bit?
                for _ in range(4):
                    MoveOneAtomTight(fixed, move, moveRange)
                    newMolecule = Atoms(buildUp, cell=boxSize)
                    newMolecule.center()
                    currentEnergy = GetEnergy(newMolecule, calc)
                    if currentEnergy < bestEnergy:
                        bestEnergy = currentEnergy
                        bestStructure = buildUp
    return bestEnergy, bestStructure


def MoveOneAtomTight(fixedAtom, atomToMove, moveRange):
    # Stop the atom drifting away into space.
    middle = fixedAtom.position
    high = moveRange * 2
    multis = random.random(3)
    directions = (multis - 0.5) * high
    atomToMove.position = middle + directions

StartEA(testC12)
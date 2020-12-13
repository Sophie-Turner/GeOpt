from Model.EAs import *


def StartEA(elementsList):
    # Set up initial values & placeholders
    # calc = SetUpVasp()
    print("EA starting with", elementsList)
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

    population = []
    bestMolecules = []
    # This will become a dataset for the potential energy plot.
    plotE = []
    plotXYZ = []

    if __name__ == 'Model.EAperAtom':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, covRads, calc) for _ in range(6)]
            for f in futures.as_completed(results):
                thisResult = f.result()
                population.append([thisResult[0], thisResult[1]])
                # Append separately from each process to prevent them being added in the wrong order by multiple process
                # attempting to write at the same time. Long-winded but done to avoid corruption and unnecessary locks.
                plotE.append(thisResult[2])
                plotXYZ.append(thisResult[3])
        population = RankByE(population, 3)
        for eachBest in population:
            newMolecule = Atoms(eachBest[0], cell=boxSize)
            newMolecule.center()
            bestMolecules.append(newMolecule)

        # newMolecule = Atoms(eachBest, pbc=false, cell=boxSize)
        # or
        # newMolecule.set_pbc((False, False, False))
        print("len(plotE):", len(plotE))
        print("len(plotXYZ):", len(plotXYZ))

    return bestMolecules, population


def Evolve(elementsList, boxSize, covRads, calc):
    allE = []
    allXYZ = []
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
        buildUpBestEnergy, buildUpBestStructure = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure, onlyH,
                                                                covRads, boxSize, calc, None, None)
    print("this best energy:", buildUpBestEnergy)

    similarity = 0
    iterations = 1
    # Wait for convergence (but don't wait for too long).
    while similarity < 2:
        # Move each atom around every atom.
        newBestEnergy, newBestStructure = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure, onlyH,
                                                        covRads, boxSize, calc, allE, allXYZ)
        print("len(allE):", len(allE))
        print("len(allXYZ):", len(allXYZ))
        if abs(newBestEnergy - buildUpBestEnergy) < 0.5:
            similarity += 1
        else:
            similarity = 0
        if newBestEnergy < buildUpBestEnergy:
            buildUpBestEnergy = newBestEnergy
            buildUpBestStructure = newBestStructure
            buildUp = newBestStructure
        iterations += 1
    return buildUp, buildUpBestEnergy, allE, allXYZ


def TestAllPlaces(buildUp, bestEnergy, bestStructure, onlyH, covRads, boxSize, calc, allE, allXYZ):
    # Move each atom around every atom.
    # This will become a dataset for the potential energy plot.
    eplotPositions = []
    eplotEnergies = []
    numSoFar = len(buildUp)
    for eachAtomToMove in range(numSoFar):
        for eachAtomFixed in range(numSoFar):
            fixed = buildUp[eachAtomFixed]
            move = buildUp[eachAtomToMove]
            # Don't let H2 form and break the molecule apart
            # but if the molecule only contains H (eg. H2) we still need to do it.
            if fixed.symbol != 'H' or move.symbol != 'H' or onlyH is True:
                moveRange = covRads[eachAtomToMove] + covRads[eachAtomFixed]
                for _ in range(6):
                    MoveOneAtomTight(fixed, move, moveRange)
                    newMolecule = Atoms(buildUp, cell=boxSize)
                    newMolecule.center()
                    currentEnergy = GetEnergy(newMolecule, calc)

                    # Now we need to maintain a dataset of all possible positions and their associated energies to
                    # show the energy plot but we only want to do this when all atoms are in the molecule.
                    if allE is not None:
                        allE.append(currentEnergy)
                        allXYZ.append(buildUp)

                    del newMolecule
                    if currentEnergy < bestEnergy:
                        bestEnergy = currentEnergy
                        bestStructure = buildUp
        # Revert to best structure at this point.
        buildUp = bestStructure
    return bestEnergy, bestStructure





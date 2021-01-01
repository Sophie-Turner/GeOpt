from Model.EAs import *


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

    population = []
    bestMolecules = []
    # This will become a dataset for the potential energy plot.
    plot = []

    if __name__ == 'Model.EAperAtom':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, covRads, calc) for _ in range(6)]
            for f in futures.as_completed(results):
                thisResult = f.result()
                population.append([thisResult[0], thisResult[1], thisResult[2]])
                # Append separately from each process to prevent them being added in the wrong order by multiple process
                # attempting to write at the same time. Long-winded but done to avoid corruption and unnecessary locks.
                plot.append(thisResult[3])
        population = RankByE(population, 3)
        for eachBest in population:
            newMolecule = Atoms(eachBest[0], cell=boxSize)
            newMolecule.center()
            bestMolecules.append(newMolecule)

        # newMolecule = Atoms(eachBest, pbc=false, cell=boxSize)
        # or
        # newMolecule.set_pbc((False, False, False))

    return bestMolecules, population, plot


def Evolve(elementsList, boxSize, covRads, calc):
    plot = []
    buildUp = []
    buildUpBestStructure = buildUp
    worstEnergy = 0
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
        buildUpBestEnergy, buildUpBestStructure, worstEnergy = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure, onlyH,
                                                                covRads, boxSize, calc, None, worstEnergy)

    similarity = 0
    iterations = 1
    # Wait for convergence (but don't wait for too long).
    while similarity < 2:
        # Move each atom around every atom.
        newBestEnergy, newBestStructure, worstEnergy = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure, onlyH,
                                                        covRads, boxSize, calc, plot, worstEnergy)
        if abs(newBestEnergy - buildUpBestEnergy) < 0.5:
            similarity += 1
        else:
            similarity = 0
        if newBestEnergy < buildUpBestEnergy:
            buildUpBestEnergy = newBestEnergy
            buildUpBestStructure = newBestStructure
            buildUp = newBestStructure
        iterations += 1
    return buildUp, buildUpBestEnergy, worstEnergy, plot


def TestAllPlaces(buildUp, bestEnergy, bestStructure, onlyH, covRads, boxSize, calc, plot, worstEnergy):
    # Move each atom around every atom.
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
                    x, y, z = MoveOneAtomTight(fixed, move, moveRange)
                    newMolecule = Atoms(buildUp, cell=boxSize)
                    newMolecule.center()
                    currentEnergy = GetEnergy(newMolecule, calc)

                    # Now we need to maintain a dataset of all possible positions and their associated energies to
                    # show the energy plot but we only want to do this when all atoms are in the molecule.
                    if plot is not None:
                        # Separate the arrays for fewer iterations and less indexing later on.
                        plot.append((currentEnergy, x, y, z, move.symbol))

                    del newMolecule
                    if currentEnergy < bestEnergy:
                        bestEnergy = currentEnergy
                        bestStructure = buildUp
                    elif currentEnergy > worstEnergy:
                        worstEnergy = currentEnergy
        # Revert to best structure at this point.
        buildUp = bestStructure
    return bestEnergy, bestStructure, worstEnergy





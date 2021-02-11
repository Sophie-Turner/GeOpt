from Model.Algos import *


def Start(elementsList, pbc, cores, numPoints):
    print('Starting Per-atom algorithm.')
    # Set up initial values & placeholders
    calc, thisPopulation, boxSize = SetUp(elementsList)
    # Get the final cell movement size
    covRads = thisPopulation.covRads

    # Move Hydrogens to the end of the list, as they can mess up the calculations at the start.
    moves = 0
    while moves < len(elementsList):
        if elementsList[0] == 'H':
            elementsList.remove('H')
            elementsList.append('H')
        else:
            break
        moves += 1

    # These will become datasets for the potential energy plots etc.
    population, bestMolecules, plot, pes, refs, energies = [], [], [], [], [], []

    if __name__ == 'Model.PerAtom':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, covRads, calc, pbc, numPoints) for _ in range(cores)]
            ProcessResults(results, population, plot, pes, refs)
        population = RankByE(population, 3)
        for eachBest in population:
            newMolecule = Atoms(eachBest[0], cell=boxSize)
            newMolecule.center()
            bestMolecules.append(newMolecule)
            energies.append([eachBest[1], eachBest[2]])

    return bestMolecules, energies, plot, pes, refs


def Evolve(elementsList, boxSize, covRads, calc, pbc, numPoints):
    plot, pes, refs, buildUp = [], [], [], []
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
        buildUpBestEnergy, buildUpBestStructure, worstEnergy = TestAllPlaces(buildUp, buildUpBestEnergy,
                                                                             buildUpBestStructure, onlyH,
                                                                             covRads, boxSize, calc, None, None, None,
                                                                             worstEnergy, pbc, numPoints)

    similarity = 0
    iterations = 1
    # Wait for convergence (but don't wait for too long).
    while similarity < 2:
        # Move each atom around every atom.
        newBestEnergy, newBestStructure, worstEnergy = TestAllPlaces(buildUp, buildUpBestEnergy, buildUpBestStructure,
                                                                     onlyH, covRads, boxSize, calc, plot, pes, refs,
                                                                     worstEnergy, pbc, numPoints)
        if abs(newBestEnergy - buildUpBestEnergy) < 0.5:
            similarity += 1
        else:
            similarity = 0
        if newBestEnergy < buildUpBestEnergy:
            buildUpBestEnergy = newBestEnergy
            buildUpBestStructure = newBestStructure
            buildUp = newBestStructure
        iterations += 1
    return buildUp, worstEnergy, buildUpBestEnergy, plot, pes, refs


def TestAllPlaces(buildUp, bestEnergy, bestStructure, onlyH, covRads, boxSize, calc, plot, pes, refs, worstEnergy, pbc, numPoints):
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
                    currentEnergy = GetEnergy(newMolecule, calc, pbc)

                    # Now we need to maintain a dataset of all possible positions and their associated energies to
                    # show the energy plot but we only want to do this when all atoms are in the molecule.
                    if plot is not None:
                        if len(plot) < numPoints:
                            # Find distances and angles for PES plot.
                            if numSoFar > 1 and eachAtomFixed != eachAtomToMove:
                                distance = newMolecule.get_distance(eachAtomToMove, eachAtomFixed)
                                if move != buildUp[-1]:
                                    otherAtom = eachAtomToMove + 1
                                    if otherAtom != eachAtomFixed:
                                        angle = newMolecule.get_angle(eachAtomToMove, eachAtomFixed, otherAtom)
                                        otherSymbol = buildUp[otherAtom].symbol
                                    else:
                                        angle = 0
                                        otherSymbol = ''
                                    reference = move.symbol + fixed.symbol + otherSymbol
                                    pes.append((distance, currentEnergy, angle, reference))
                                    if reference not in refs:
                                        refs.append(reference)
                            # Separate the arrays for fewer iterations and less indexing later on.
                            plot.append((currentEnergy, x, y, z, move.symbol))
                    # Clean up to avoid leaving molecules lying around all over the place.
                    del newMolecule
                    if currentEnergy < bestEnergy:
                        bestEnergy = currentEnergy
                        bestStructure = buildUp
                    elif currentEnergy > worstEnergy:
                        worstEnergy = currentEnergy
        # Revert to best structure at this point.
        buildUp = bestStructure
    return bestEnergy, bestStructure, worstEnergy





from Model.Algos import *
from time import time


def StartEA(elementsList, pbc, popSize, cores, numPoints, mutDist, mutSize, permute, cross):
    print('Starting Many-molecule EA with a population size of', popSize)
    # Set up and initialise our template molecule to start with.
    calc, thisPopulation, boxSize = SetUp(elementsList, mutSize)
    print('box size', boxSize)
    bestMolecules, energies, plot, pes, refs, finalists = [], [], [], [], [], []
    firstCoordinates = thisPopulation.initPositions
    atomObjectList = thisPopulation.initAtomsObject
    numAtoms = len(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()
    parentEnergy = GetEnergy(parentMolecule, calc, pbc)

    if numAtoms == 1:
        # If there's only 1 atom we can skip all this...
        pass

    # Start with a reasonably spaced estimate to use as an energy comparison.
    population = [[parentMolecule, firstCoordinates, parentEnergy]]
    lastBestEnergy = parentEnergy

    # Start with some permutations.
    for i in range(popSize):
        # Make some permutations.
        MakeNewMolecule(elementsList, firstCoordinates, lastBestEnergy, None, boxSize, population,
                        False, False, True, None, None, calc, pbc, numPoints)
    # Find the best permutation.
    if cores > 1:
        population = RankByE(population, 2)

    # Use multiprocessing to quickly compare results.
    if __name__ == 'Model.EAmanyMolecules':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, population, cores, calc, pbc, popSize, numPoints,
                                       mutDist, cross, permute)for _ in range(cores)]
            ProcessResults(results, finalists, plot, pes, refs)
    if cores > 2:
        finalists = RankByE(finalists, 3)
    elif cores == 2:
        finalists = RankByE(finalists, 2)
    for eachMolecule in finalists:
        bestMolecules.append(eachMolecule[0])
        energies.append([eachMolecule[1], eachMolecule[2]])

    print('EA has finished')
    return bestMolecules, energies, plot, pes, refs


def Evolve(elementsList, boxSize, population, cores, calc, pbc, popSize, numPoints, mutDist, cross, permute):
    print('Evolving the populations')
    # These will be private to each thread.
    plot, pes = [], []
    worstEnergy = 0
    # Ranges of random atom movements.
    width = boxSize[0]
    changeSizes = [width/30, width/20, width/16, width/12, width/8, width/4]
    # Set the mutation distributions.
    mutMove, mutFill = mutDist+1, mutDist+3
    # See how many iterations it takes.
    iterations = 0
    similarity = 0

    bestCoordinates = population[0][1]
    lastBestEnergy = population[0][2]
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/startLoop.png", population[0][0],
         rotation='10x,30y,0z')
    population.pop(0)

    # End if the best energy doesn't change much for several consecutive iterations.
    while similarity < 5 and iterations < 20:

        # New child molecules.
        for i in range(popSize):
            # Make some with random mutations.
            _, _ = MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, changeSizes[randint(6)], boxSize,
                                   population, mutMove, False, permute, plot, pes, calc, pbc, numPoints)
        # Selection.
        if cores > 1:
            population = RankByE(population, 2)
        bestCoordinates = population[0][1]

        # Create some random molecules.
        for i in range(popSize):
            # Introduce some random strangers.
            outRefs, outNewEnergy = MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, None, boxSize,
                                                    population, mutFill, cross, permute, plot, pes, calc, pbc, numPoints)
            if outRefs is not None:
                refs, newEnergy = outRefs, outNewEnergy
                if newEnergy > worstEnergy:
                    worstEnergy = newEnergy

        # Selection.
        if cores > 1:
            population = RankByE(population, 2)

        # Update the stopping criterion.
        newBestEnergy = population[0][2]

        if abs(lastBestEnergy - newBestEnergy) < 0.25:
            similarity = similarity + 1
        else:
            similarity = 0
        iterations = iterations + 1

    for member in population:
        bestMolecule = member[0]

    print(iterations, 'iterations')

    return bestMolecule, worstEnergy, newBestEnergy, plot, pes, refs

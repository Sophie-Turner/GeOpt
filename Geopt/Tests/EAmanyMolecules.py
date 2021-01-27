from Model.Algos import *
from time import time


def StartEA(elementsList):
    startTime = time()
    # Set up and initialise our template molecule to start with.
    # calc = SetUpVasp()
    calc = EMT()
    thisPopulation = Population(elementsList)
    bestMolecules, returnPop, plot, pes, refs = [], [], [], [], []
    boxSize = thisPopulation.boxSize
    firstCoordinates = thisPopulation.initPositions
    atomObjectList = thisPopulation.initAtomsObject
    numAtoms = len(elementsList)

    parentMolecule = Atoms(atomObjectList, cell=boxSize)
    parentMolecule.center()
    parentEnergy = GetEnergy(parentMolecule, calc)

    if numAtoms == 1:
        # If there's only 1 atom we can skip all this...
        pass

    # Start with a reasonably spaced estimate to use as an energy comparison.
    population = [[parentMolecule, firstCoordinates, parentEnergy]]
    lastBestEnergy = population[0][2]

    # Start with some permutations.
    for i in range(len(elementsList)):
        # Make some permutations.
        _, _, _, _ = MakeNewMolecule(elementsList, firstCoordinates, lastBestEnergy, None, boxSize, population,
                                     False, False, True, calc)

    # Find the best permutation.
    population = RankByE(population, 1)

    if __name__ == '__main__':
        with futures.ProcessPoolExecutor() as executor:
            results = [executor.submit(Evolve, elementsList, boxSize, population, calc) for _ in range(6)]
            for f in futures.as_completed(results):
                thisResult = f.result()
                bestMolecules.append(thisResult[0])
                returnPop.append(thisResult[1])
                plot.append(thisResult[2])
                pes.append(thisResult[3])
                refs.append(thisResult[4])

    endTime = time()
    print('Time taken = {} seconds'.format(round(endTime-startTime)))

    return bestMolecules, returnPop, plot, pes, refs


def Evolve(elementsList, boxSize, population, calc):
    # This will be a private dataset to each thread.
    returnPop = []
    worstEnergy = 0
    # Ranges of random atom movements.
    width = boxSize[0]
    changeSizes = [width/30, width/20, width/16, width/12, width/8, width/4]
    # See how many iterations it takes.
    iterations = 0
    similarity = 0

    bestCoordinates = population[0][1]
    lastBestEnergy = population[0][2]
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/startLoop.png", population[0][0],
         rotation='10x,30y,0z')
    population.pop(0)

    # End if the best energy doesn't change much for several consecutive iterations.
    while similarity < 5 and iterations < 100:
        # New child molecules.
        for i in range(len(elementsList)):
            # Make some with random mutations.
            _, _, _, _ = MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, changeSizes[1], boxSize,
                                         population, 2, False, False, calc)

        # Selection.
        population = RankByE(population, 1)
        bestCoordinates = population[0][1]

        # Create some random molecules.
        for i in range(len(elementsList)):
            # Introduce some random strangers.
            plot, pes, refs, newEnergy = MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, None, boxSize,
                                                         population, 3, False, False, calc)
            if newEnergy > worstEnergy:
                worstEnergy = newEnergy

        # Selection.
        population = RankByE(population, 1)

        # Update the stopping criterion.
        newBestEnergy = population[0][2]

        if abs(lastBestEnergy - newBestEnergy) < 0.25:
            similarity = similarity + 1
        else:
            similarity = 0
        iterations = iterations + 1

    for member in population:
        bestMolecule = member[0]
        # Put this outside the loop to reduce the number of []s.
        theseAtoms = []
        for atom in bestMolecule:
            theseAtoms.append(atom)
        returnPop.append([theseAtoms, worstEnergy, newBestEnergy])

    return bestMolecule, returnPop, plot, pes, refs


elementsList = ['H', 'H', 'O']
one, two, three, four, five = StartEA(elementsList)
print('bestMolecules:', one)
print('Number of best molecules:', len(one))
print('population:', two)
print('plot:', three)
print('pes:', four)
print('refs', five)


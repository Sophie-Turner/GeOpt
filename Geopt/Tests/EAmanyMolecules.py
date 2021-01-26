from Model.Algos import *
from time import time

global bestMolecules, returnPop, plot, pes, refs

def StartEA(elementsList):
    # Set up and initialise our template molecule to start with.
    #calc = SetUpVasp()
    startTime = time()
    calc = EMT()
    thisPopulation = Population(elementsList)
    global bestMolecules, returnPop, plot, pes, refs
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
        MakeNewMolecule(elementsList, firstCoordinates, lastBestEnergy, None, boxSize, population, False, False,
                        True, plot, pes, refs, calc)

    # Find the best permutation.
    population = RankByE(population, 1)

    Evolve(elementsList, boxSize, population, calc)

    endTime = time()
    print('Time taken = {} seconds'.format(round(endTime-startTime)))

    return bestMolecules, returnPop, plot, pes, refs


def Evolve(elementsList, boxSize, population, calc):
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
            MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, changeSizes[1], boxSize, population, 2,
                            False, False, plot, pes, refs, calc)

        # Selection.
        population = RankByE(population, 1)
        bestCoordinates = population[0][1]

        # Create some random molecules.
        for i in range(len(elementsList)):
            # Introduce some random strangers.
            MakeNewMolecule(elementsList, bestCoordinates, lastBestEnergy, None, boxSize, population, 3, False,
                            False, plot, pes, refs, calc)

        # Selection.
        population = RankByE(population, 3)

        # Update the stopping criterion.
        newBestEnergy = population[0][2]

        if abs(lastBestEnergy - newBestEnergy) < 0.25:
            similarity = similarity + 1
        else:
            similarity = 0
        iterations = iterations + 1

    print("Iterations performed: ", iterations)
    print("Population size: ", len(population))
    print('Population:', population)

    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/optimised.png", population[0][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/optimised2.png", population[1][0],
          rotation='10x,30y,0z')
    write("C:/Users/pipin/Documents/fyp/SophieCOMP3000/Geopt/Images/optimised3.png", population[2][0],
          rotation='10x,30y,0z')

    for member in population:
        molecule = member[0]
        theseAtoms = []
        for atom in molecule:
            theseAtoms.append(atom)
        bestMolecules.append(molecule)
        returnPop.append(theseAtoms)


elementsList = ['H', 'H', 'O']
one, two, three, four, five = StartEA(elementsList)
print('bestMolecules:', one)
print('population:', two)
print('plot:', three)
print('pes:', four)
print('refs', five)


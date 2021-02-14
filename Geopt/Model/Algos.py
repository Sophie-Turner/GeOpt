from ase import Atoms, Atom
from ase.io import write
from ase.calculators.vasp import Vasp
from Model.EmtCalculator import EMT
from Model.Populations import Population
from numpy import random
from concurrent import futures


def SetUp(elementsList, mutSize):
    print('setting up calculator and cell')
    # calc = SetUpVasp()
    calc = EMT()
    thisPopulation = Population(elementsList, mutSize)
    boxSize = thisPopulation.boxSize
    return calc, thisPopulation, boxSize


def SetUpVasp():
    calc = Vasp(
        xc='pbe',  # exchange-correlation functional
        nbands=6,  # number of bands
        encut=350,  # planewave cutoff
        ismear=1,  # Methfessel-Paxton smearing
        sigma=0.01,  # very small smearing factor for a molecule
    )
    return calc


def ProcessResults(results, population, plot, pes, refs):
    for f in futures.as_completed(results):
        thisResult = f.result()
        population.append([thisResult[0], thisResult[1], thisResult[2]])
        plot.append(thisResult[3])
        pes.append(thisResult[4])
        refs.append(thisResult[5])
    print('processing results from multiprocessing')


def MakeNewMolecule(elementsList, inCoordinates, bestE, changeSize, boxSize, population, mutate, cross, permute, plot,
                    pes, calc, pbc, numPoints):
    if cross is True:
        inCoordinates = Crossover(population)
    if permute is True:
        inCoordinates = random.permutation(inCoordinates)
    childCoordinates, atomsObject = PrepareChild(elementsList, inCoordinates)
    if mutate == 1:
        MoveAtomsUniform(atomsObject, changeSize)
    elif mutate == 2:
        MoveAtomsGauss(atomsObject, 0, changeSize)
    elif mutate == 3:
        FillCellUniform(atomsObject, boxSize)
    elif mutate == 4:
        FillCellGauss(atomsObject, boxSize)
    elif mutate == 5:
        MoveHydrogen(atomsObject, boxSize)
    childMolecule = GenerateChild(atomsObject, boxSize)
    childEnergy = GetEnergy(childMolecule, calc, pbc)
    # Don't keep any extremely terrible structures.
    if childEnergy >= bestE * 100:
        print('Energy too high. Recursing')
        MakeNewMolecule(elementsList, inCoordinates, bestE, changeSize, boxSize, population, mutate, cross, permute,
                        plot, pes, calc, pbc, numPoints)
    population.append([childMolecule, childCoordinates, childEnergy])
    if plot is not None:
        if len(plot) < numPoints:
            refs = []
            numAtoms = len(childMolecule)
            for each in range(numAtoms):
                eachAtom = childMolecule[each]
                coords = eachAtom.position
                x, y, z = coords[0], coords[1], coords[2]
                plot.append((childEnergy, x, y, z, eachAtom.symbol))
                for other in range(numAtoms-1):
                    if each != other:
                        otherAtom = childMolecule[other]
                        distance = childMolecule.get_distance(each, other)
                        if numAtoms > 2 and each != other+1:
                            nextAtom = childMolecule[other + 1]
                            angle = childMolecule.get_angle(each, other, other+1)
                            ref = (eachAtom.symbol + otherAtom.symbol + nextAtom.symbol)
                        else:
                            angle = 0
                            ref = (eachAtom.symbol + otherAtom.symbol)
                        pes.append((distance, childEnergy, angle, ref))
                        if ref not in refs:
                            refs.append(ref)
            return refs, childEnergy
        else:
            return None, childEnergy


def PrepareChild(elementsList, startCoordinates):
    childAtomsObject = []
    # Need to copy the parent's co-ordinates so we don't alter the parent molecule at the same time.
    childCoordinates = startCoordinates[:]
    for i in range(len(elementsList)):
        childAtomsObject.append(Atom(elementsList[i], childCoordinates[i]))
    return childCoordinates, childAtomsObject


def GenerateChild(childAtomsObject, boxSize):
    child = Atoms(childAtomsObject, cell=boxSize)
    # Translate atoms to the centre of the unit cell. It's OK if the atoms stick out of their cell.
    child.center()
    return child


def GetEnergy(molecule, calc, pbc):
    # Set up the ase force calculator for finding energies.
    molecule.set_pbc(pbc)
    molecule.calc = calc
    energy = molecule.get_potential_energy()
    return energy


def MoveOneAtomTight(fixedAtom, atomToMove, moveRange):
    # Stop the atom drifting away into space.
    middle = fixedAtom.position
    high = moveRange * 2
    multis = random.random(3)
    directions = (multis - 0.5) * high
    atomToMove.position = middle + directions
    coords = atomToMove.position
    x, y, z = coords[0], coords[1], coords[2]
    return x, y, z


def MoveOneAtomGauss(fixedAtom, atomToMove, sigma):
    mean = fixedAtom.position
    atomToMove.position = (random.normal(mean, sigma, 3))


def MoveAtomsUniform(itsAtoms, changeSize):
    for eachAtom in itsAtoms:
        eachAtom.position += ((random.randint(-changeSize-1, changeSize+1, 3)) / 100)


def MoveAtomsGauss(itsAtoms, mean, sigma):
    for eachAtom in itsAtoms:
        eachAtom.position += (random.normal(mean, sigma, 3))


def FillCellUniform(itsAtoms, boxSize):
    for eachAtom in itsAtoms:
        eachAtom.position = (random.random(3)*boxSize/2)


def FillCellGauss(itsAtoms, boxSize):
    centre = boxSize[0]/2
    for eachAtom in itsAtoms:
        eachAtom.position = (random.normal(centre, centre/3, 3))


def PermuteAtoms(itsAtoms, inCoordinates):
    newCoordinates = random.permutation(inCoordinates)
    for eachAtom in itsAtoms:
        eachAtom.position = newCoordinates[eachAtom]


def MoveHydrogen(itsAtoms, boxSize):
    # Hydrogens tend to get lost so push them all onto the Carbon atoms.
    change = boxSize[0]/12
    moveTo = []
    carbons = 0
    for eachAtom in itsAtoms:
        if eachAtom.symbol == 'C':
            moveTo.append(eachAtom.position)
        elif eachAtom.symbol == 'H':
            if len(moveTo) > 0:
                try:
                    eachAtom.position = moveTo[carbons] + (random.normal(0, change, 3))
                    carbons += 1
                except:
                    carbons = 0
                    eachAtom.position = moveTo[carbons] + (random.normal(0, change, 3))
            else:
                eachAtom.position += (random.normal(0, change * 2, 3))


def RankByE(population, numToKeep):
    # I originally planned to use a quicksort but decided that since the list was small
    # there was no need to make it complicated.
    rankedPopulation = []
    bestMolecule = None
    while len(rankedPopulation) < numToKeep:
        bestEnergy = 1000
        for eachMember in population:
            # Access the last item of population because different algorithms have different types
            # of population and put low energy at the end.
            if eachMember[-1] < bestEnergy:
                bestEnergy = eachMember[-1]
                bestMolecule = eachMember
        rankedPopulation.append(bestMolecule)
        population.remove(bestMolecule)
    return rankedPopulation


def Crossover(population):
    mumCoordinates = population[0][1]
    howMany = len(mumCoordinates)
    mumOrDad = random.randint(0, 2, howMany)
    childCoordinates = []
    for i in range(howMany):
        whichParent = mumOrDad[i]
        thisPoint = population[whichParent][1][i]
        childCoordinates.append(thisPoint)
    return childCoordinates


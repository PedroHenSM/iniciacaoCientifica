#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
@author pedromuniz
"""
# init-hook="from pylint.config import find_pylintrc; import os, sys; sys.path.append(os.path.dirname(find_pylintrc()))"

import sys
sys.path.append("../functions")
sys.path.append("../functions/utils")
sys.path.append("../functions/cec2020/boundConstrained")
sys.path.append("../functions/eureka")
import utils
import cec2020BoundConstrained
import eureka
import simpleFunctions as Functions
import numpy as np
import timeit
import operator as op # For sorting population
from copy import deepcopy


# Constants
DEB = "DEB"
APM = "APM"



class Individual(object):
  def __init__(self, n, objectiveFunction, g=None, h=None, violations=None, violationSum=None, fitness=None):
    self.n = n
    self.objectiveFunction = objectiveFunction
    self.g = g
    self.h = h
    self.violations = violations
    self.violationSum = violationSum
    self.fitness = fitness

  def printIndividual(self, boolObjFunc, constraintHandling, boolN):
    if boolObjFunc:
      print("{}".format(self.objectiveFunction[0]), end=" ")
    if constraintHandling is not None:
      if constraintHandling == DEB:
        print("{}".format(self.violationSum), end=" ")
      elif constraintHandling == APM:
        print("{}".format(self.fitness), end=" ")
    if boolN:
      print(*self.n, sep=" ")

  # Makes print(individual) a string, not a reference (memory adress)
  def __repr__(self):
    return str(self.__dict__) + "\n"

class Population(object):
  def __init__(self, nSize, popSize, function, objFunctionSize, lowerBound, upperBound, gSize, hSize):
    strFunction = str(function)
    self.individuals = []
    for _ in range(popSize):
      n = []
      objFunc = [99999 for i in range(objFunctionSize)]
      g = [99999 for i in range(gSize)] if gSize is not None else None
      h = [99999 for i in range(hSize)] if hSize is not None else None
      violations = [99999 for i in range(gSize + hSize)] if gSize is not None else None
      violationSum = 99999 if gSize is not None else None
      for _ in range(nSize):
        if strFunction[0] == "1": # Truss problems
          n.append(np.random.uniform(lowerBound, upperBound))
        elif strFunction[0] == "3": # cec 2020 bound constrained
          n.append(np.random.uniform(-100, 100))
        else:
          sys.exit("Function not defined.")
      self.individuals.append(Individual(n, objFunc, g, h, violations, violationSum))

  def printPopulation(self, boolObjFunc, constraintHandling, boolN):
    for individual in self.individuals:
      if boolObjFunc:
        print(individual.objectiveFunction)
      if constraintHandling:
        if constraintHandling == DEB:
          print(individual.violationSum)
        if constraintHandling == APM:
          print(individual.fitness)
      if boolN:
        print(individual.n)
    print()

  def printBest(self, boolObjFunc, constraintHandling, boolN):
    best = self.bestIndividual(constraintHandling)
    best.printIndividual(boolObjFunc, constraintHandling, boolN)

  def evaluate(self, function, fe, truss, case, discreteSet):
    strFunction = str(function)
    objFuncSize = len(self.individuals[0].objectiveFunction)
    nSize = len(self.individuals[0].n)
    for individual in self.individuals:
      fe = fe + 1
      if strFunction[0] == "1": # Trusses
        # Values array/list is the array/list with the objectiveFunction and constraints of type g
        valuesArraySize = truss.getNumberObjectives() + truss.getNumberConstraints() # the size will be objFunction (1) + gSize
        xArray = utils.new_doubleArray(truss.getDimension()) # creates an array
        valuesArray = utils.new_doubleArray(valuesArraySize) # the size will be objFunct(1) + gSize
        # Transfers values from a python list to a C++ array
        populateArray(xArray, individual.n, 0, truss.getDimension(), case, discreteSet)
        valuesList = individual.objectiveFunction + individual.g
        populateArray(valuesArray, valuesList, 0, valuesArraySize, case, discreteSet)
        # Evaluate population
        truss.evaluation(xArray, valuesArray)
        # Transfers values from a C++ array to a python list
        # populateList(individual.n, xArray, 0, truss.getDimension()) #TODO Verificar se pode ficar comentado
        individual.objectiveFunction[0] = utils.doubleArray_getitem(valuesArray, 0)
        populateList(individual.g, valuesArray, 1, valuesArraySize)
        # Cleans mess
        utils.delete_doubleArray(xArray)
        utils.delete_doubleArray(valuesArray)
      elif strFunction[0] == "2":
        sys.exit("Not defined")
      elif strFunction[0] == "3": # Cec 2020 bound constrained
        # Gets all numbers except the first
        func = int(strFunction[1:])
        # Creates an array with nSize dimension for project variables and objectiveFunction
        xArray = utils.new_doubleArray(nSize)
        objFuncArray = utils.new_doubleArray(objFuncSize)
        # Transfers values from a python list to a C++ array
        populateArray(xArray, individual.n, 0, nSize, case, discreteSet)
        populateArray(objFuncArray, individual.objectiveFunction, 0, objFuncSize, case, discreteSet)
        # Evaluate population
        cec2020BoundConstrained.cec20_test_func(xArray, objFuncArray, nSize, 1, func)
        # Transfers values from a C++ array to a python list
        populateList(individual.n, xArray, 0, nSize)
        populateList(individual.objectiveFunction, objFuncArray, 0, objFuncSize)
        # Cleans mess
        utils.delete_doubleArray(xArray)
        utils.delete_doubleArray(objFuncArray)
      else:
        sys.exit("Function not defined.")

    # Returns functions evaluations
    return fe

  def deSelect(self, offsprings, generatedOffspring, constraintHandling):
    idx = 0
    for i in range(len(self.individuals)): # Parents size
      bestIdx = idx
      # walks through every n offsprings of each parent
      while idx < generatedOffspring * (i + 1):
        # Picks the best offspring
        if constraintHandling is None: # Bound constrained problems
          if offsprings.individuals[idx].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
            bestIdx = idx
        elif constraintHandling == DEB:
          if offsprings.individuals[idx].violationSum < offsprings.individuals[bestIdx].violationSum:
            bestIdx = idx
          elif offsprings.individuals[idx].violationSum == offsprings.individuals[bestIdx].violationSum:
            if offsprings.individuals[idx].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
              bestIdx = idx
        elif constraintHandling == APM:
          if offsprings.individuals[idx].fitness < offsprings.individuals[bestIdx].fitness:
            bestIdx = idx
          elif offsprings.individuals[idx].fitness == offsprings.individuals[bestIdx].fitness:
            if offsprings.individuals[idx].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
              bestIdx = idx
        else:
          sys.exit("Constraint handling method not defined.")
        idx+=1

      # Best offsprings better than parent, he becomes the parent
      if constraintHandling is None: # Bound constrained problems
        if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
          self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
      elif constraintHandling == DEB:
        if offsprings.individuals[bestIdx].violationSum < self.individuals[i].violationSum:
          self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
        elif offsprings.individuals[bestIdx].violationSum == self.individuals[i].violationSum:
          if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
            self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
      elif constraintHandling == APM:
        if offsprings.individuals[bestIdx].fitness < self.individuals[i].fitness:
          self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
        elif offsprings.individuals[bestIdx].fitness == self.individuals[i].fitness:
          if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
            self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])

  def checkBounds(self, function, lowerBound, upperBound):
    strFunction = str(function)
    nMin = nMax = 0
    if strFunction[0] == "1": # Truss problems
      nMin = lowerBound
      nMax = upperBound
    elif strFunction[0] == "3": # Cec 2020 bound constrained functions
      nMin = -100
      nMax = 100
    else:
      sys.exit("Function not defined.")
    for individual in self.individuals:
      for i in range(len(individual.n)):
        if individual.n[i] > nMax:
          individual.n[i] = nMax
        elif individual.n[i] < nMin:
          individual.n[i] = nMin

  def sort(self, offsprings, constraintHandling):
    if constraintHandling is None: # Bound constrained problem
      self.individuals.sort(key=op.attrgetter("objectiveFunction"))
      if offsprings is not None:
        offsprings.individuals.sort(key=op.attrgetter("objectiveFunction"))
    elif constraintHandling == DEB:
      self.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
      if offsprings is not None:
        offsprings.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
    elif constraintHandling == APM:
      self.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))
      if offsprings is not None:
        offsprings.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))
    else:
      sys.exit("Constraint handling not defined.")

  def bestIndividual(self, constraintHandling):
    # Assume the first individual is the best one
    best = deepcopy(self.individuals[0])
    for individual in self.individuals:
      if constraintHandling is None: # Bound constrained problems
        if individual.objectiveFunction[0] < best.objectiveFunction[0]:
          best = deepcopy(individual)
      elif constraintHandling == DEB:
        if individual.violationSum < best.violationSum:
          best = deepcopy(individual)
        elif individual.violationSum == best.violationSum:
          if individual.objectiveFunction[0] < best.objectiveFunction[0]:
            best = deepcopy(individual)
      elif constraintHandling == APM:
        bothViolates = individualViolates(individual, constraintHandling) and individualViolates(best, constraintHandling) # Both violates
        neitherViolates = not individualViolates(individual, constraintHandling) and not individualViolates(best, constraintHandling) # Neither violates
        if bothViolates or neitherViolates:
          if individual.fitness < best.fitness:
            best = deepcopy(individual)
          elif individual.fitness == best.fitness:
            if individual.objectiveFunction[0] < best.objectiveFunction[0]:
              best = deepcopy(individual)
        elif individualViolates(best, constraintHandling): # Best violates
          best = deepcopy(individual)
      else:
        sys.exit("Constraint handling method not defined.")

      # elif constraintHandling == 2: # APM
      # 	if self.individuals[i].fitness < best.fitness:
      # 		best = deepcopy(self.individuals[i])
      # 	elif self.individuals[i].fitness == best.fitness:
      # 		if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
      # 			best = deepcopy(self.individuals[i])

    return best

  def bestIndividualDeprecated(self, constraintHandling):
    # Assume the first individual is the best one
    best = deepcopy(self.individuals[0])
    for i in range(1, len(self.individuals)):
      if constraintHandling == DEB:
        if self.individuals[i].violationSum < best.violationSum:
          best = deepcopy(self.individuals[i])
        elif self.individuals[i].violationSum == best.violationSum:
          if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
            best = deepcopy(self.individuals[i])
      elif constraintHandling == APM:
        bothViolates = individualViolates(self.individuals[i], constraintHandling) and individualViolates(best, constraintHandling) # Both violates
        neitherViolates = not individualViolates(self.individuals[i], constraintHandling) and not individualViolates(best, constraintHandling) # Neither violates
        if bothViolates or neitherViolates:
          if self.individuals[i].fitness < best.fitness:
            best = deepcopy(self.individuals[i])
          elif self.individuals[i].fitness == best.fitness:
            if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
              best = deepcopy(self.individuals[i])
        elif individualViolates(best, constraintHandling): # Best violates
          best = deepcopy(self.individuals[i])

      # elif constraintHandling == 2: # APM
      # 	if self.individuals[i].fitness < best.fitness:
      # 		best = deepcopy(self.individuals[i])
      # 	elif self.individuals[i].fitness == best.fitness:
      # 		if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
      # 			best = deepcopy(self.individuals[i])
      else:
        if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
          best = deepcopy(self.individuals[i])
    return best

  def hallOfFame(self, hof, constraintHandling):
    currentBest = self.bestIndividual(constraintHandling)
    if hof is None:
      hof = deepcopy(currentBest)
    if constraintHandling is None: # Bound constrained problems
      if currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
        hof = deepcopy(currentBest)
    elif constraintHandling == DEB: # Deb
      if currentBest.violationSum < hof.violationSum:
        hof = deepcopy(currentBest)
      elif currentBest.violationSum == hof.violationSum:
        if currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
          hof = deepcopy(currentBest)
    elif constraintHandling == APM: # APM
      bothViolates = individualViolates(currentBest, constraintHandling) and individualViolates(hof, constraintHandling) # Both violates
      neitherViolates = not individualViolates(currentBest, constraintHandling) and not individualViolates(hof, constraintHandling) # Neither violates
      if bothViolates or neitherViolates:
        if currentBest.fitness < hof.fitness:
          hof = deepcopy(currentBest)
        elif currentBest.fitness == hof.fitness:
          if currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
            hof = deepcopy(currentBest)
      elif individualViolates(hof, constraintHandling): # Hof violates
        hof = deepcopy(currentBest)
    else:
      sys.exit("Constraint handling method not defined.")

      # if hof is None or currentBest.fitness < hof.fitness:
      # 	hof = deepcopy(currentBest)
      # elif hof is None or currentBest.fitness == hof.fitness:
      # 	if hof is None or currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
      # 		hof = deepcopy(currentBest)

    return hof

  def uniteConstraints(self, constraintHandling):
    gSize = len(self.individuals[0].g)
    hSize = len(self.individuals[0].h)
    constraintsSize = gSize + hSize
    for individual in self.individuals:
      idxG = idxH = 0
      for i in range(constraintsSize):
        if i<gSize:
          individual.violations[i] = individual.g[idxG]
          idxG+=1
        else:
          if constraintHandling == DEB:
            individual.violations[i] = individual.h[idxH]
          elif constraintHandling == APM:
            individual.violations[i] = np.abs(individual.h[idxH]) - 0.0001
          idxH+=1
      if constraintHandling == DEB:
        # Only sums positives values
        individual.violationSum = np.sum(value for value in individual.violations if value > 0)

  def calculatePenaltyCoefficients(self, numberOfConstraints, penaltyCoefficients, averageObjectiveFunctionValues):
    popSize = len(self.individuals)
    sumObjectiveFunction = 0
    # foreach candidate solution
    for individual in self.individuals:
      sumObjectiveFunction += individual.objectiveFunction[0]
    # the absolute value of sumObjectiveFunction
    np.abs(sumObjectiveFunction)
    # the average of the objective function values
    averageObjectiveFunctionValues = sumObjectiveFunction / popSize
    # the denominator of the equation of the penalty coefficients
    denominator = 0
    # the sum of the constraint violation values
    # these values are recorded to be used in the next situation

    sumViolation = []
    for i in range(numberOfConstraints):
      sumViolation.append(0)
      for individual in self.individuals:
        if individual.violations[i] > 0:
          sumViolation[i] += individual.violations[i]
      denominator += sumViolation[i] *  sumViolation[i]
    
    # the penalty coefficients are calculated
    for i in range(numberOfConstraints):
      if denominator == 0:
        penaltyCoefficients[i] = 0
      else:
        penaltyCoefficients[i] = (sumObjectiveFunction / denominator) * sumViolation[i]

    return averageObjectiveFunctionValues

  def calculateAllFitness(self, numberOfConstraints, penaltyCoefficients, averageObjectiveFunctionValues):
    for individual in self.individuals:
      # indicates if the candidate solution is feasible
      infeasible = False
      # the penalty value
      penalty = 0
      for i in range(numberOfConstraints):
        if individual.violations[i] > 0:
          # the candidate solution is infeasible if some constraint is violated
          infeasible = True
          # the penalty value is updated
          penalty += penaltyCoefficients[i] * individual.violations[i]
      # fitness is the sum of the objective function and penalty values
      # the candidate solution is infeasible and just the objective function value,
      # otherwise
      if infeasible:
        if individual.objectiveFunction[0] > averageObjectiveFunctionValues:
          individual.fitness = individual.objectiveFunction[0] + penalty
        else:
          individual.fitness = averageObjectiveFunctionValues + penalty
      else:
        individual.fitness = individual.objectiveFunction[0]

  def handleConstraints(self, population, constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc):
    self.uniteConstraints(constraintHandling)
    avgObjFunc = self.calculatePenaltyCoefficients(constraintsSize, penaltyCoefficients, avgObjFunc)
    self.calculateAllFitness(constraintsSize, penaltyCoefficients, avgObjFunc)
    # TODO Verificar se é necessário recalcular o fitness dos pais após avgObjFunc modificar (creio que seja)
    if population is not None:
      population.calculateAllFitness(constraintsSize, penaltyCoefficients, avgObjFunc)
    return avgObjFunc

  def moveArzToPop(self, list2d):
    for index, individual in enumerate(self.individuals):
      individual.n = deepcopy(list2d[index])

  def selectMuIndividuals(self, mu):
    population = []
    for i in range(mu):
      population.append(self.individuals[i].n)
    return population

  def deGeneratePopulation(self, offsprings, generatedOffspring, CR, F):
    parentsSize = len(self.individuals)
    nSize = len(self.individuals[0].n)
    offspringIdx = 0
    for i in range(parentsSize):
      for _ in range(generatedOffspring):
        chosenOnesIdxs = populationPick(i, parentsSize)
        R = np.random.randint(0, nSize)  # Random index
        for j in range(nSize):
          Ri = np.random.rand()  # Generates random number between (0,1)
          if Ri < CR or j == R:
            offsprings.individuals[offspringIdx].n[j] = self.individuals[chosenOnesIdxs[0]].n[j] + F * (self.individuals[chosenOnesIdxs[1]].n[j] - self.individuals[chosenOnesIdxs[2]].n[j])
          else:
            offsprings.individuals[offspringIdx].n[j] = self.individuals[i].n[j]
        offspringIdx = offspringIdx + 1

  def cmaGeneratePopulation(self, parentsSize, centroid, sigma, BD):
    # Generates new array of project variables
    nSize = len(self.individuals[0].n)
    arz = np.random.standard_normal((parentsSize, nSize))
    arz = centroid + sigma * np.dot(arz, BD.T)
    l2d = arz.tolist() # Transforms matrix in a list of lists
    self.moveArzToPop(l2d)

  def cmaUpdateCovarianceMatrix(self, parentsSize, mu, centroid, sigma, weights, mueff, cc, cs, ccov1, ccovmu, damps, pc, ps, B, diagD, C, update_count, chiN,BD):
    nSize = len(self.individuals[0].n)
    populationList2d = self.selectMuIndividuals(mu)
    old_centroid = centroid
    centroid = np.dot(weights, populationList2d) # Recombination

    c_diff = centroid - old_centroid

    # Cumulation : update evolution paths
    ps = (1 - cs) * ps \
        + np.sqrt(cs * (2 - cs) * mueff) / sigma \
        * np.dot(B, (1. / diagD) *
                    np.dot(B.T, c_diff))

    hsig = float((np.linalg.norm(ps) /
                  np.sqrt(1. - (1. - cs) ** (2. * (update_count + 1.))) / chiN <
                  (1.4 + 2. / (nSize + 1.))))

    update_count += 1

    pc = (1 - cc) * pc + hsig \
        * np.sqrt(cc * (2 - cc) * mueff) / sigma \
        * c_diff

    # Update covariance matrix
    artmp = populationList2d - old_centroid
    C = (1 - ccov1 - ccovmu + (1 - hsig) *
              ccov1 * cc * (2 - cc)) * C \
        + ccov1 * np.outer(pc, pc) \
        + ccovmu * np.dot((weights * artmp.T), artmp) \
        / sigma ** 2

    # Adapt step-size sigma
    sigma *= np.exp((np.linalg.norm(ps) / chiN - 1.) *
                            cs / damps)

    diagD, B = np.linalg.eigh(C)
    indx = np.argsort(diagD)

    cond = diagD[indx[-1]] / diagD[indx[0]]

    diagD = diagD[indx] ** 0.5
    B = B[:, indx]
    BD = B * diagD
    return centroid, sigma, pc, ps, B, diagD, C, update_count, BD

# Print(population) is shown as string and not memory reference
  def __repr__(self):
    return str(self.__dict__) + "\n"

def printInitialPopulationInfo(algorithm, constraintHandling, function, seed, parentsSize, offspringsSize, maxFe, CR, F, rweights, case):
  feasibiliyMeasure = "-"
  if constraintHandling is None:
    pass
  elif constraintHandling == DEB:
    feasibiliyMeasure = "ViolationSum"
  elif constraintHandling == APM:
    feasibiliyMeasure = "Fitness"
  else:
    sys.exit("Constraint handling method not defined.")
  print("Algorithm: {}".format(algorithm))
  print("Constrainth handling method: {}".format(constraintHandling))
  print("Function: {}".format(function))
  print("Case: {}".format(case))
  print("Seed: {}".format(seed))
  print("Parents size: {}".format(parentsSize))
  print("Offsprings size: {}".format(offspringsSize))
  print("maxFe: {}".format(maxFe))
  print("CR: {}".format(CR))
  print("F: {}".format(F))
  print("rweights: {}".format(rweights))
  print("*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*")
  print("ngen ObjectiveFunction {} ProjectVariables".format(feasibiliyMeasure))

def printFinalPopulationInfo(status, population, hof, constraintHandling):
  print("*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*")
  print("Status: {}".format(status))
  print("Last individual")
  population.printBest(True, constraintHandling, True)
  print("Hall of fame")
  hof.printIndividual(True, constraintHandling, True)

def initializeTruss(function):
  if function == 110:  # Truss 10 bars
    truss = eureka.F101Truss10Bar()
  elif function == 125:  # Truss 25 bars
    truss = eureka.F103Truss25Bar()
  elif function == 160:  # Truss 60 bars
    truss = eureka.F105Truss60Bar()
  elif function == 172:  # Truss 72 bars
    truss = eureka.F107Truss72Bar()
  elif function == 1942:  # Truss 942 bars
    truss = eureka.F109Truss942Bar()
  else:
    sys.exit("Function not defined.")
  bounds = utils.new_doubleddArray(truss.getDimension())
  bounds = utils.castToDouble(truss.getBounds())
  lowerBound = utils.doubleddArray_getitem(bounds, 0, 0)
  upperBound = utils.doubleddArray_getitem(bounds, 0, 1)
  return truss, lowerBound, upperBound

def getDiscreteCaseList(function):
  if function == 110:
    chosenSet = [1.62, 1.80, 1.99, 2.13, 2.38, 2.62, 2.63, 2.88, 2.93, 3.09,
      3.13, 3.38, 3.47, 3.55, 3.63, 3.84, 3.87, 3.88, 4.18, 4.22, 4.49, 4.59,
      4.80, 4.97, 5.12, 5.74, 7.22, 7.97, 11.50, 13.50, 13.90, 14.20, 15.50,
      16.00, 16.90, 18.80, 19.90, 22.00, 22.90, 26.50, 30.00, 33.50]
  elif function == 125:
    chosenSet = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.1, 3.2]
  elif function == 160:
    chosenSet = [0.5, 0.6, 0.7, 0.8, 0.9, 1. , 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. , 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3. , 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4. , 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.]
  elif function == 172:
    chosenSet = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. , 1.1, 1.2, 1.3,
        1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. , 2.1, 2.2, 2.3, 2.4, 2.5]
  elif function == 1942:
    chosenSet = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
        12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22.,
        23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33.,
        34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44.,
        45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55.,
        56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66.,
        67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77.,
        78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88.,
        89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
       100., 101., 102., 103., 104., 105., 106., 107., 108., 109., 110.,
       111., 112., 113., 114., 115., 116., 117., 118., 119., 120., 121.,
       122., 123., 124., 125., 126., 127., 128., 129., 130., 131., 132.,
       133., 134., 135., 136., 137., 138., 139., 140., 141., 142., 143.,
       144., 145., 146., 147., 148., 149., 150., 151., 152., 153., 154.,
       155., 156., 157., 158., 159., 160., 161., 162., 163., 164., 165.,
       166., 167., 168., 169., 170., 171., 172., 173., 174., 175., 176.,
       177., 178., 179., 180., 181., 182., 183., 184., 185., 186., 187.,
       188., 189., 190., 191., 192., 193., 194., 195., 196., 197., 198.,
       199., 200.]
  return chosenSet

# Find nearestvalue in list (considering that list is sorted)
def findNearest(array,value):
  idx = np.searchsorted(array, value, side="left")
  if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
    return array[idx-1]
  else:
    return array[idx]

def find_nearest(array, value):
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return array[idx]

def initiliazeHandleConstraintsParams(function):
  gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize = initializeConstraints(function)
  penaltyCoefficients = [None for i in range(constraintsSize)]
  avgObjFunc = -1
  return gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc

def individualViolates(individual, constraintHandling):
  if constraintHandling == DEB:
    if individual.violationSum == 0:
      return False
  elif constraintHandling == APM:
    if individual.objectiveFunction[0] == individual.fitness:
      return False
  else:
    sys.exit("Constriaint handling method not defined.")
  return True

def constraintsInitParams (function, constraintHandling):
  if constraintHandling is not None:
    gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc = initiliazeHandleConstraintsParams(function)
  else:
    sys.exit("Constraint handling not defined")

  return gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc

def cmaInitParams(nSize, parentsSize, centroid, sigma, mu, rweights):
  # User defined params
  if parentsSize is None:
    parentsSize = int(4 + 3 * np.log(nSize)) # Population size
  if centroid is None:
    # centroid = np.array([5.0]*nSize) # objective variables initial points
    centroid = np.random.randn(nSize) # objective variables initial points | gaussian distribution between [0,1]
  if sigma is None:
    # sigma = 5.0 # coordinate wise standard deviation(step-size)
    sigma = 0.5 # coordinate wise standard deviation(step-size)
  if mu is None:
    mu = int(parentsSize / 2) # number of parents selected (selected search points in the population)
  if rweights is None:
    rweights = "linear"

  # Initiliaze dynamic (internal) strategy parameters and constants
  pc = np.zeros(nSize) # evolution paths for C
  ps = np.zeros(nSize) # evolutions paths for sigma
  chiN = np.sqrt(nSize) * (1 - 1. / (4. * nSize) + 1. / (21. * nSize ** 2)) # expectation of ||N(0,I)|| == norm(randn(N,1))
  C = np.identity(nSize) # covariance matrix
  diagD, B = np.linalg.eigh(C) # 
  indx = np.argsort(diagD) # return the indexes that would sort an array
  diagD = diagD[indx] ** 0.5 # diagonal matrix D defines the scaling
  B = B[:, indx] # B defines de coordinate system
  BD = B * diagD
  cond = diagD[indx[-1]] / diagD[indx[0]] # divide o valor mais alto do pelo mais baixo
  update_count = 0 # B and D update at feval == 0

  # Strategy parameter setting: Selection
  if rweights == "superlinear":
    weights = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1))
  elif rweights == "linear":
    weights = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1)) # muXone recombination weights
  elif rweights == "equal":
    weights = np.ones(mu)
  weights /= np.sum(weights) # normalize recombination weights array
  mueff = 1. / np.sum(weights **2) # vriance-effective size of mu

  # Strategy parameter setting: Adaptation
  cc = 4. / (nSize + 4.) # time constant for cumulation for C
  cs = (mueff + 2.) / (nSize + mueff + 3.) # t-const for cumulation for sigma control
  ccov1 = 2. / ((nSize + 1.3) ** 2 + mueff) # learning rate for rank one update of C
  ccovmu = 2. * (mueff - 2. + 1. / mueff) / ((nSize + 2.) ** 2 + mueff) # and for rank-mu update
  ccovmu = min(1 - ccov1, ccovmu)
  damps = 1. + 2. * max(0, np.sqrt((mueff - 1.) / (nSize + 1.)) - 1.) + cs # # damping for sigma

  return parentsSize, mu, centroid, sigma, pc, ps, chiN, C, diagD, B, BD, update_count, weights, mueff, cc, cs, ccov1, ccovmu, damps

def deInitParams(function, parentsSize, offspringsSize):
  # Define DE parameters
  CR = 0.9
  F = 0.6
  generatedOffspring = int(offspringsSize / parentsSize)
  return CR, F, generatedOffspring

# Python function to pass values of a python list to a C++ (swig) array
def populateArray(a, l, startIdx, size, case, discreteSet):
  for i in range(size - startIdx):
    # Sets the arr[startIdx] with value from list[i]
    if case == "continuous":
      utils.doubleArray_setitem(a, startIdx, l[i])
    elif case == "discrete":
      # utils.doubleArray_setitem(a, startIdx, findNearest(discreteSet, l[i]))
      utils.doubleArray_setitem(a, startIdx, float(find_nearest(discreteSet, l[i])))
    else:
      sys.exit("Case not discrete or continuous.")
    startIdx+=1

# Python function to pass values of a C++ array (swig) to a python list
def populateList(l, a, startIdx, size):
  for i in range(size - startIdx):
    # Sets the list[i] with the value from arr[startIdx]
    l[i] = utils.doubleArray_getitem(a, startIdx)
    startIdx+=1

def defineMaxEval(function, nSize):
  strFunction = str(function)
  maxFe = None
  if strFunction[0] == "2": # Trusses
    sys.exit("Not implemented")
  elif strFunction[0] == "3": # cec2020 bound constrained
    if(nSize == 5):
      maxFe = 50000
    elif (nSize == 10):
      maxFe = 1000000
    elif (nSize == 15):
      maxFe = 3000000
    elif (nSize == 20):
      maxFe = 10000000
    else:
      sys.exit("Dimension not defined for cec2020 bound constrained functions.")
  else:
    sys.exit("Function not defined.")

  return maxFe

# Generate random indexes that won't repeat on a generation of new offsprings
def populationPick(parentIdx, parentsSize):
  chosenOnes = []
  while len(chosenOnes) < 3:
    idx = np.random.randint(0, parentsSize)
    if idx != parentIdx:
      if idx not in chosenOnes:
        chosenOnes.append(idx)
  return chosenOnes

def initializeConstraints(function):
  strFunction = str(function)
  truss = lowerBound = upperBound = None
  if strFunction[0] == "1": # Trusses
    truss, lowerBound, upperBound = initializeTruss(function)
    nSize = truss.getDimension()
    g, h, constraintsSize = truss.getNumberConstraints(), 0, truss.getNumberConstraints() + 0
  else:
    sys.exit("Constraints not defined for function.")
  return g, h, constraintsSize, truss, lowerBound, upperBound, nSize

# Differential evolution
def DE(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case):
  np.random.seed(seed)
  strFunction = str(function)
  feval = 0
  hof = None
  status = "Initializing"

  if strFunction[0] == "3": # cec2020 bound constrained
    maxFe = defineMaxEval(function, nSize)

  if strFunction[0] != "1": # Bound constrained problems
    constraintHandling = None

  # Initialize differential evolution parameters
  CR, F, generatedOffspring = deInitParams(function, parentsSize, offspringsSize)

  # Initialize truss constraints, if necessary
  if constraintHandling:
    gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc = constraintsInitParams(function, constraintHandling)
    discreteSet = getDiscreteCaseList(function)
  else:
    lowerBound = upperBound = gSize = hSize = truss = None

  # Generate initial population and evaluate
  parents = Population(nSize, parentsSize, function, 1, lowerBound, upperBound, gSize, hSize)
  offsprings = Population(nSize, offspringsSize, function, 1, lowerBound, upperBound, gSize, hSize)
  feval = parents.evaluate(function, feval, truss, case, discreteSet)

  # Constraints handling
  if constraintHandling:
    avgObjFunc = parents.handleConstraints(None, constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)

  # DE+DEB10keval: T10: 5071.401883658249 | T25: 484.0522202840531 | T60: 428.33801862127643 | T72: 391.41071625298974 
  # DE+APM10keval: T10: 5067.045222753923 | T25: 484.065205988603 | T60: 425.8179547877767 | T72: 383.963083800768
  printInitialPopulationInfo("Differential Evolution", constraintHandling, function, seed, parentsSize, offspringsSize, maxFe, CR, F, "-", case)

  if feval > maxFe:
    sys.exit("Maximum number of function evaluations too low.")
  while feval < maxFe:
    status = "Executing"
    print(feval, end=" ")
    # Generate new population
    parents.deGeneratePopulation(offsprings, generatedOffspring, CR, F)

    # Check bounds and evaluate offsprings
    offsprings.checkBounds(function, lowerBound, upperBound)
    feval = offsprings.evaluate(function, feval, truss, case, discreteSet)
    # Handling constraints, if necessary
    if constraintHandling:
      avgObjFunc = offsprings.handleConstraints(parents, constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)

    # Selects best offsprings and move them into parents
    parents.deSelect(offsprings, generatedOffspring, constraintHandling)

    # Gets hall of fame individual
    hof = parents.hallOfFame(hof, constraintHandling)
    
    # Prints best individual of current generation
    parents.printBest(True, constraintHandling, True)
  status = "Finished"
  printFinalPopulationInfo(status, parents, hof, constraintHandling)

# CMA ES
def CMAES(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case):
  np.random.seed(seed)
  strFunction = str(function)
  feval = 0
  hof = None
  status="Initializing"

  if strFunction[0] == "3": # cec2020 bound constrained
    maxFe = defineMaxEval(function, nSize)

  if strFunction[0] != "1": # Bound constrained problems
    constraintHandling = None


  # User defined params (if set to None, uses default)
  parentsSize = centroid = sigma = mu = rweights = None
  # rweights = "equal"

  if constraintHandling:
    gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc = constraintsInitParams(function, constraintHandling)
    discreteSet = getDiscreteCaseList(function)
  else:
    lowerBound = upperBound = gSize = hSize = truss = None

  # Initialize all cmaes params
  parentsSize, mu, centroid, sigma, pc, ps, chiN, C, diagD, B, BD, update_count, weights, mueff, cc, cs, ccov1, ccovmu, damps = cmaInitParams(nSize, parentsSize, centroid, sigma, mu, rweights)

  # Generate initial population
  parents = Population(nSize, parentsSize, function, 1, lowerBound, upperBound, gSize, hSize)


  # CMAESrweights=equal+DEB: T10: 5061.970157299801 | T25: 484.05229151214365 | T60: 311.6048214519769 | T72: 380.0802294286406  
  # CMAES+APM: T10: 5062 | T25: 484 | T60: 313 | T72: 380
  printInitialPopulationInfo("CMA-ES", constraintHandling, function, seed, parentsSize, mu, maxFe, "-", "-", rweights, case)

  if feval > maxFe:
    sys.exit("Maximum number of function evaluations too low.")
  while feval < maxFe:
    status = "Executing"
    print(feval, end=" ")
    # Generate new population
    parents.cmaGeneratePopulation(parentsSize, centroid, sigma, BD)
    # Check bounds and evaluate
    parents.checkBounds(function, lowerBound, upperBound)
    feval = parents.evaluate(function, feval, truss, case, discreteSet)


    # Handling constraints, if necessary
    if constraintHandling:
      avgObjFunc = parents.handleConstraints(None, constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)
  
    # Sorts population
    parents.sort(None, constraintHandling)

    # Gets hall of fame individual
    hof = parents.hallOfFame(hof, constraintHandling)

    # Update covariance matrix strategy from the population
    centroid, sigma, pc, ps, B, diagD, C, update_count, BD = parents.cmaUpdateCovarianceMatrix(
      parentsSize, mu, centroid, sigma, weights, mueff, cc, cs, ccov1, ccovmu, damps, pc, 
      ps, B, diagD, C, update_count, chiN, BD
    )

    # Prints best individual of current generation
    parents.printBest(True, constraintHandling, True)
  status = "Finished"
  # Prints final info
  printFinalPopulationInfo(status, parents, hof, constraintHandling)

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
import engineeringProblems
import numpy as np
import timeit
import operator as op # For sorting population
import warnings
from collections import deque # Double ended queue
warnings.filterwarnings("error", category=RuntimeWarning)
from copy import deepcopy
# np.seterr(all='raise')

# Constants
DEB = "DEB"
APM = "APM"
TOLFUN = 10**-6
TOLX = 10**-12

AlGORITHM_PARAMS = {
  "function": None,
  "parentsSize": None,
  "offspringsSize": None,
  "seed": None,
  "maxFe": None,
  "constraintHandling": None,
  "case": None,
  "rweights": None
}



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
    self.cmaRestart = False
    for _ in range(popSize):
      n = []
      objFunc = [99999 for i in range(objFunctionSize)]
      g = [99999 for i in range(gSize)] if gSize is not None else None
      h = [99999 for i in range(hSize)] if hSize is not None else None
      violations = [99999 for i in range(gSize + hSize)] if gSize is not None else None
      violationSum = 99999 if gSize is not None else None
      for i in range(nSize):
        if strFunction[0] == "1": # Truss problems
          n.append(np.random.uniform(lowerBound, upperBound))
        elif strFunction[0] == "2": # Mechanical Engineering Problems
          if function == 21: # The Tension/Compression Spring Design
            if i == 0:
              n.append(np.random.uniform(2, 15))
            elif i == 1:
              n.append(np.random.uniform(0.25, 1.3))
            elif i == 2:
              n.append(np.random.uniform(0.05, 2))
            else:
              sys.exit("Design variable should not exist.")
          elif function == 22: # The Speed Reducer design
            if i == 0:
              n.append(np.random.uniform(2.6, 3.6))
            elif i == 1:
              n.append(np.random.uniform(0.7, 0.8))
            elif i == 2:
              n.append(np.round(np.random.uniform(17, 28)))
            elif i == 3:
              n.append(np.random.uniform(7.3, 8.3))
            elif i == 4:
              n.append(np.random.uniform(7.8, 8.3))
            elif i == 5:
              n.append(np.random.uniform(2.9, 3.9))
            elif i == 6:
              n.append(np.random.uniform(5, 5.5))
            else:
              sys.exit("Design variable should not exist.")
          elif function == 222: # The Speed Reducer design (ObjFunc 2994)
            if i == 0:
              n.append(np.random.uniform(2.6, 3.6))
            elif i == 1:
              n.append(np.random.uniform(0.7, 0.8))
            elif i == 2:
              n.append(np.round(np.random.uniform(17, 28)))
            elif i == 3:
              n.append(np.random.uniform(7.3, 8.3))
            elif i == 4:
              n.append(np.random.uniform(7.3, 8.3))
            elif i == 5:
              n.append(np.random.uniform(2.9, 3.9))
            elif i == 6:
              n.append(np.random.uniform(5, 5.5))
            else:
              sys.exit("Design variable should not exist.")
          elif function == 23: # The Welded Beam design
            if i == 0:
              n.append(np.random.uniform(0.125, 10))
            elif i == 1:
              n.append(np.random.uniform(0.1, 10))
            elif i == 2:
              n.append(np.random.uniform(0.1, 10))
            elif i == 3:
              n.append(np.random.uniform(0.1, 10))
            else:
              sys.exit("Design variable should not exist.")
          elif function == 233: # The Welded Beam design (ObjFunc 1.7)
            if i == 0:
              n.append(np.random.uniform(0.1, 2))
            elif i == 1:
              n.append(np.random.uniform(0.1, 10))
            elif i == 2:
              n.append(np.random.uniform(0.1, 10))
            elif i == 3:
              n.append(np.random.uniform(0.1, 2))
            else:
              sys.exit("Design variable should not exist.")
          elif function == 24: # The Pressure Vessel design
            if i == 0:
              n.append(np.random.uniform(0.0625, 5))
            elif i == 1:
              n.append(np.random.uniform(0.0625, 5))
            elif i == 2:
              n.append(np.random.uniform(10, 200))
            elif i == 3:
              n.append(np.random.uniform(10, 200))
            else:
              sys.exit("Design variable should not exist.")
          elif function == 25: # The Cantilever Beam design
            # n[0] -> h1 | n[1] -> b1 | n[2] -> h2 | n[3] -> b2 (...)
            if i == 0: # b1
              n.append(np.round(np.random.uniform(1, 5)))
            elif i == 1: # h1
              n.append(np.round(np.random.uniform(30, 65)))
            elif i == 2: # b2
              n.append(np.random.uniform(2.4, 3.1))
            elif i == 3: # h2
              n.append(np.random.uniform(45, 60))
            elif i == 4: # b3
              n.append(np.random.uniform(2.4, 3.1))
            elif i == 5: # h3
              n.append(np.random.uniform(45, 60))
            elif i == 6: # b4
              n.append(np.random.uniform(1, 5))
            elif i == 7: # h4
              n.append(np.random.uniform(30, 65))
            elif i == 8: # b5
              n.append(np.random.uniform(1, 5))
            elif i == 9: # h5
              n.append(np.random.uniform(30, 65))
            else:
              sys.exit("Design variable should not exist.")
            # elif i == 10:
            #   sys.exit("Verify this last constraint. Not implemented.")
          else:
            sys.exit("Function not defined.")

        elif strFunction[0] == "3": # cec 2020 bound constrained
          n.append(np.random.uniform(-100, 100))
        else:
          sys.exit("Function not defined.")
      # 60 bar truss, best solution from Erica
      # n = [1.16186, 2.05396, 0.55622, 1.56836, 1.42066, 0.56823, 1.92526, 1.90424, 1.11336, 1.49336, 1.49575, 0.80284, 2.06031, 1.30629, 1.00834, 0.54040, 0.50147, 0.94529, 1.14428, 1.19860, 1.15119, 1.08666, 0.50457, 1.02406, 1.26823]
      # 60 bar truss, best solution from Silva
      # n = [1.1628, 2.0741, 0.5000, 1.7550, 1.6937, 0.5831, 1.9866, 1.8817, 1.0210, 1.8557, 1.7946, 0.5001, 2.0715, 1.2532, 1.0677, 0.7102, 0.7795, 1.0434, 1.1508, 1.1600, 0.9856, 1.0803, 0.6282, 1.0245, 1.2642]
      # 60 bar truss, best solution from Bernardino
      # n = [1.1726, 2.0621, 0.5100, 1.7654, 1.6658, 0.5720, 1.9062, 1.9484, 1.0510, 1.7466, 1.6457, 0.5153, 2.0797, 1.2616, 1.1398, 0.6888, 0.7888, 1.0401, 1.1605, 1.1499, 1.0097, 1.0680, 0.8199, 1.0693, 1.2743]
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
        # print("valuesArraySize: {}".format(valuesArraySize))
        # print("xArraySize: {}".format(truss.getDimension()))
        # sys.exit('ok');
        # # TEST
        # nErica = [1.16186, 2.05396, 0.55622, 1.56836, 1.42066, 0.56823, 1.92526, 1.90424, 1.11336, 1.49336, 1.49575, 0.80284, 2.06031, 1.30629, 1.00834, 0.54040, 0.50147, 0.94529, 1.14428, 1.19860, 1.15119, 1.08666, 0.50457, 1.02406, 1.26823]
        # USE_T25D
        # nErica = [10.904000, 16.802600, 15.436700, 5.304100, 4.921300, 13.856200, 16.115300, 7.447200]
        # nErica = [0.1, 0.3, 3.4, 0.1, 2.1, 1.0, 0.5, 3.4 ]
        # USE_T10C
        # nErica = [12.658400, 19.451960, 17.878820, 6.208860, 5.767980, 16.058520, 18.660380, 8.677120, 8.473380, 10.410580]
        # USE_T60C
        # nErica = [2.192000, 3.107300, 2.895350, 1.323050, 1.263650, 2.650100, 3.000650, 1.655600, 1.628150, 1.889150, 2.327900, 3.731000, 4.477100, 1.138550, 4.968500, 4.868600, 2.436350, 3.672500, 2.619500, 0.883400, 2.613200, 0.783050, 4.551800, 4.419050, 0.888800 ]
        # populateArray(xArray, nErica, 0, truss.getDimension(), case, discreteSet)
        # #  END TEST
        # Transfers values from a python list to a C++ array
        populateArray(xArray, individual.n, 0, truss.getDimension(), case, discreteSet)
        valuesList = individual.objectiveFunction + individual.g
        populateArray(valuesArray, valuesList, 0, valuesArraySize, case, discreteSet)
        # Evaluate population
        truss.evaluation(xArray, valuesArray)
        # Transfers values from a C++ array to a python list
        # populateList(individual.n, xArray, 0, truss.getDimension()) # Verificar se pode ficar comentado
        individual.objectiveFunction[0] = utils.doubleArray_getitem(valuesArray, 0)
        populateList(individual.g, valuesArray, 1, valuesArraySize)
        # print("obj func: {}".format(utils.doubleArray_getitem(valuesArray, 0)))
        # for i in range(truss.getNumberConstraints()):
        #   print("{} ".format(individual.g[i]))
        # sys.exit("saindo no evaluate")
        # Cleans mess
        utils.delete_doubleArray(xArray)
        utils.delete_doubleArray(valuesArray)
      elif strFunction[0] == "2":
        individual.objectiveFunction[0], individual.g, individual.h = engineeringProblems.executeFunction(function, individual.n, individual.objectiveFunction, individual.g, individual.h)
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
    nMinList = nMaxList = None
    nMin = nMax = 0
    if strFunction[0] == "1": # Truss problems
      nMin = lowerBound
      nMax = upperBound
    elif strFunction[0] == "2": # Engineering problems
      if function == 21: # The Tension/Compression Spring Design
        nMinList = [2, 0.25, 0.05]
        nMaxList = [15, 1.3, 2]
      elif function == 22: # The Speed Reducer design
        nMinList = [2.6, 0.7, 17, 7.3, 7.8, 2.9, 5]
        nMaxList = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5]
      elif function == 222: # The Speed Reducer design (ObjFunc 2994)
        nMinList = [2.6, 0.7, 17, 7.3, 7.3, 2.9, 5]
        nMaxList = [3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5]
      elif function == 23: # The Welded Beam design
        nMinList = [0.125, 0.1, 0.1, 0.1]
        nMaxList = [10, 10, 10, 10]
      elif function == 233: # The Welded Beam design (ObjFunc 1.7)
        nMinList = [0.1, 0.1, 0.1, 0.1]
        nMaxList = [2, 10, 10, 2]
      elif function == 24: # The Pressure Vessel design
        nMinList = [0.0625, 0.0625, 10, 10]
        nMaxList = [5, 5, 200, 200]
      elif function == 25: # The Cantilever Beam design
        nMinList = [1, 30, 2.4, 45, 2.4, 45, 1, 30, 1, 30]
        nMaxList = [5, 65, 3.1, 60, 3.1, 60, 5, 65, 5, 65]
      else:
        sys.exit("Function not defined.")
    elif strFunction[0] == "3": # Cec 2020 bound constrained functions
      nMin = -100
      nMax = 100
    else:
      sys.exit("Function not defined.")
    for individual in self.individuals:
      for i in range(len(individual.n)):
        if nMaxList is not None and nMinList is not None: # Each dimension has a different bound
          nMin = nMinList[i]
          nMax = nMaxList[i]
        # Verify if each design variable respects the bounds
        if individual.n[i] < nMin:
          individual.n[i] = nMin
        elif individual.n[i] > nMax:
          individual.n[i] = nMax

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
        bothViolates = isInfactible(individual, constraintHandling) and isInfactible(best, constraintHandling) # Both violates
        neitherViolates = not isInfactible(individual, constraintHandling) and not isInfactible(best, constraintHandling) # Neither violates
        if bothViolates or neitherViolates:
          if individual.fitness < best.fitness:
            best = deepcopy(individual)
          elif individual.fitness == best.fitness:
            if individual.objectiveFunction[0] < best.objectiveFunction[0]:
              best = deepcopy(individual)
        elif isInfactible(best, constraintHandling): # Best violates
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
      bothViolates = isInfactible(currentBest, constraintHandling) and isInfactible(hof, constraintHandling) # Both violates
      neitherViolates = not isInfactible(currentBest, constraintHandling) and not isInfactible(hof, constraintHandling) # Neither violates
      if bothViolates or neitherViolates:
        if currentBest.fitness < hof.fitness:
          hof = deepcopy(currentBest)
        elif currentBest.fitness == hof.fitness:
          if currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
            hof = deepcopy(currentBest)
      elif isInfactible(hof, constraintHandling): # Hof violates
        hof = deepcopy(currentBest)
    else:
      sys.exit("Constraint handling method not defined.")

      # if hof is None or currentBest.fitness < hof.fitness:
      # 	hof = deepcopy(currentBest)
      # elif hof is None or currentBest.fitness == hof.fitness:
      # 	if hof is None or currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
      # 		hof = deepcopy(currentBest)

    return hof


  def hallOfFameFixed(self, hof, bestLastFactible, constraintHandling):
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
      bothViolates = isInfactible(currentBest, constraintHandling) and isInfactible(hof, constraintHandling) # Both violates
      neitherViolates = not isInfactible(currentBest, constraintHandling) and not isInfactible(hof, constraintHandling) # Neither violates
      if bothViolates or neitherViolates:
        if currentBest.fitness < hof.fitness:
          hof = deepcopy(currentBest)
        elif currentBest.fitness == hof.fitness:
          if currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
            hof = deepcopy(currentBest)
      elif isInfactible(hof, constraintHandling): # Hof violates
        hof = deepcopy(currentBest)
    else:
      sys.exit("Constraint handling method not defined.")

      # if hof is None or currentBest.fitness < hof.fitness:
      # 	hof = deepcopy(currentBest)
      # elif hof is None or currentBest.fitness == hof.fitness:
      # 	if hof is None or currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
      # 		hof = deepcopy(currentBest)

    return hof

  def lastFactibleIndividual(self, lastFactible, constraintHandling):
    currentBest = self.bestIndividual(constraintHandling)
    # Current best is factible
    if not isInfactible(currentBest, constraintHandling):
      # First factible individual
      # if lastFactible is None:
      lastFactible = deepcopy(currentBest)
      # In all cases, just compare objective function, because both individuals (currentBest and bestLastFactible) are already factible
      # if constraintHandling is None or constraintHandling == DEB or constraintHandling == APM:
        # if currentBest.objectiveFunction[0] < bestLastFactible.objectiveFunction[0]:
          # bestLastFactible = deepcopy(currentBest)
      # else:
      #   sys.exit("Constraint handling method not defined.")
    return lastFactible
    

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

  def cmaUpdateCovarianceMatrix(self, parentsSize, mu, centroid, sigma, weights, mueff, cc, cs, ccov1, ccovmu, damps, pc, ps, B, diagD, C, update_count, chiN, BD):
    # nSize = len(self.individuals[0].n)
    # populationList2d = self.selectMuIndividuals(mu)
    # old_centroid = centroid
    # centroid = np.dot(weights, populationList2d) # Recombination

    # c_diff = centroid - old_centroid

    # # Cumulation : update evolution paths
    # ps = (1 - cs) * ps \
    #     + np.sqrt(cs * (2 - cs) * mueff) / sigma \
    #     * np.dot(B, (1. / diagD) *
    #                 np.dot(B.T, c_diff))

    # hsig = float((np.linalg.norm(ps) /
    #               np.sqrt(1. - (1. - cs) ** (2. * (update_count + 1.))) / chiN <
    #               (1.4 + 2. / (nSize + 1.))))

    # update_count += 1

    # pc = (1 - cc) * pc + hsig \
    #     * np.sqrt(cc * (2 - cc) * mueff) / sigma \
    #     * c_diff

    # # Update covariance matrix
    # artmp = populationList2d - old_centroid
    # C = (1 - ccov1 - ccovmu + (1 - hsig) *
    #           ccov1 * cc * (2 - cc)) * C \
    #     + ccov1 * np.outer(pc, pc) \
    #     + ccovmu * np.dot((weights * artmp.T), artmp) \
    #     / sigma ** 2

    # # Adapt step-size sigma
    # sigma *= np.exp((np.linalg.norm(ps) / chiN - 1.) *
    #                         cs / damps)

    # diagD, B = np.linalg.eigh(C)
    # indx = np.argsort(diagD)

    # cond = diagD[indx[-1]] / diagD[indx[0]]

    # diagD = diagD[indx] ** 0.5
    # B = B[:, indx]
    # BD = B * diagD
    # return centroid, sigma, pc, ps, B, diagD, C, update_count, BD
    try:
      # print("sellf.printsbest INICIO TRY: {}".format(self.printBest(True, AlGORITHM_PARAMS["constraintHandling"], True)))
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
      # print("diagD: {}".format(type(diagD)))
      # print("diagD[indx]: {}".format(diagD[indx]))
      # diagD = diagD[indx] ** 0.5
      diagD = np.power(diagD[indx], 0.5)
      # print("não CONTINUA")
      B = B[:, indx]
      BD = B * diagD
      # print("sellf.printsbest FINAL TRY: {}".format(self.printBest(True, AlGORITHM_PARAMS["constraintHandling"], True)))
      return centroid, sigma, pc, ps, B, diagD, C, update_count, BD, cond
    except (RuntimeWarning, RuntimeError, TypeError, np.linalg.LinAlgError):
    # except:
      print("Covariance matrix did not converge. Restarting algorithm")
      self.cmaRestart = True

  # Returns true if entire population is factible
  def isAllFactible(self, constraintHandling):
    for individual in self.individuals:
      if isInfactible(individual, constraintHandling):
        return False
    return True

  # Returns true if entire population is infactible
  def isAllInfactible(self, constraintHandling):
    for individual in self.individuals:
      if not isInfactible(individual, constraintHandling):
        return False
    return True

  def restartCriterias(self, constraintHandling, restartCriterias, populationDeque, noEffectAxisIdx, centroid, sigma, cond, pc, B, diagD, C):
    # restartAlgorithm = False
    # Check conditions for tolFun
    # Stop if range of the best objective function values of the last 10 + ⌈30n/λ⌉ generations is zero (equalfunvalhist) or the range of these functions values and all functions values of the recent generation is below TolFun
    # Last 10 + ⌈30n/λ⌉ objective functions on array, verify tolFun criteria
    if len(populationDeque) == populationDeque.maxlen:
      isPopulationDequeFactible = False
      isAllPopulationFactible = self.isAllFactible(constraintHandling)
      isAllPopulationInfactible = self.isAllInfactible(constraintHandling)

      # Create list for current generation
      currentGen = []
      for individual in self.individuals:
        currentGen.append(individual)
        
      # Bound constrained problems or factible individuals, check only objective function
      if constraintHandling is None or isAllPopulationFactible:
        maxCurrentGen = max(currentGen, key=lambda x: min(x.objectiveFunction))
        minCurrentGen = min(currentGen, key=lambda x: min(x.objectiveFunction))
        maxPopulationDeque = max(populationDeque, key=lambda x: min(x.objectiveFunction))
        minPopulationDeque = min(populationDeque, key=lambda x: min(x.objectiveFunction))
        print("CurrentGen Range tolFun (AllFactible): {} - {} = {}".format(maxCurrentGen.objectiveFunction[0], minCurrentGen.objectiveFunction[0], maxCurrentGen.objectiveFunction[0] - minCurrentGen.objectiveFunction[0]))
        print("PopulationDeque Range tolFun (AllFactible): {} - {} = {}".format(maxPopulationDeque.objectiveFunction[0], minPopulationDeque.objectiveFunction[0], maxPopulationDeque.objectiveFunction[0] - minPopulationDeque.objectiveFunction[0]))
        # Returns true if range from all objective functions of current generation is below tolFun
        isCurrentRangeBelowTolFun = True if maxCurrentGen.objectiveFunction[0] - minCurrentGen.objectiveFunction[0] < TOLFUN else False
        # Returns true if range from last 10 + ⌈30n/λ⌉ is below tolFun
        isRangeBelowTolFun = True if maxPopulationDeque.objectiveFunction[0] - minPopulationDeque.objectiveFunction[0] < TOLFUN else False
        # If both conditions are true, algorithm should be restarted
        restartCriterias['tolFun'] = True if isRangeBelowTolFun and isCurrentRangeBelowTolFun else False
      # Entire population is infactible, check which penalty method is being used
      elif isAllPopulationInfactible:
        if constraintHandling == DEB:
          maxCurrentGen = max(currentGen, key=lambda x: x.violationSum)
          minCurrentGen = min(currentGen, key=lambda x: x.violationSum)
          maxPopulationDeque = max(populationDeque, key=lambda x: x.violationSum)
          minPopulationDeque = min(populationDeque, key=lambda x: x.violationSum)
          print("CurrentGen Range tolFun DEB: {} - {} = {}".format(maxCurrentGen.violationSum, minCurrentGen.violationSum, maxCurrentGen.violationSum - minCurrentGen.violationSum))
          print("PopulationDeque Range tolFun DEB: {} - {} = {}".format(maxPopulationDeque.violationSum, minPopulationDeque.violationSum, maxPopulationDeque.violationSum - minPopulationDeque.violationSum))
          # Returns true if range from all objective functions of current generation is below tolFun`
          isCurrentRangeBelowTolFun = True if maxCurrentGen.violationSum - minCurrentGen.violationSum < TOLFUN else False
          # Returns true if range from last 10 + ⌈30n/λ⌉ is below tolFun
          isRangeBelowTolFun = True if maxPopulationDeque.violationSum - minPopulationDeque.violationSum < TOLFUN else False
          # If both conditions are true, algorithm should be restarted
          restartCriterias['tolFun'] = True if isRangeBelowTolFun and isCurrentRangeBelowTolFun else False
        elif constraintHandling == APM:
          maxCurrentGen = max(currentGen, key=lambda x: x.fitness)
          minCurrentGen = min(currentGen, key=lambda x: x.fitness)
          maxPopulationDeque = max(populationDeque, key=lambda x: x.fitness)
          minPopulationDeque = min(populationDeque, key=lambda x: x.fitness)
          print("CurrentGen Range tolFun APM: {} - {} = {}".format(maxCurrentGen.fitness, minCurrentGen.fitness, maxCurrentGen.fitness - minCurrentGen.fitness))
          print("PopulationDeque Range tolFun APM: {} - {} = {}".format(maxPopulationDeque.fitness, minPopulationDeque.fitness, maxPopulationDeque.fitness - minPopulationDeque.fitness))
          # Returns true if range from all objective functions of current generation is below tolFun`
          isCurrentRangeBelowTolFun = True if maxCurrentGen.fitness - minCurrentGen.fitness < TOLFUN else False
          # Returns true if range from last 10 + ⌈30n/λ⌉ is below tolFun
          isRangeBelowTolFun = True if maxPopulationDeque.fitness - minPopulationDeque.fitness < TOLFUN else False
          # If both conditions are true, algorithm should be restarted
          restartCriterias['tolFun'] = True if isRangeBelowTolFun and isCurrentRangeBelowTolFun else False
      else:
        print("Population alternating between factibles and infactibles individuals!")

    # Check tolX condition
    # Stop if the standard deviation of the normal distribution is smaller than in all coordinates and sigma*pc is amller than TolX
    # All components of pc and sqrt(diag(C)) (standard deviation of the normal distribution) are smaller than the threshold
    if all(pc < TOLX) and all(np.sqrt(np.diag(C)) < TOLX):
      # All components of pc and sqrt(diag(C)) are smaller than the threshold
      restartCriterias["tolX"] = True

    # Check noEffectAxis
    # Stop if adding 0.1-standard deviation vector in any principal axis direction of C doesnt change m
    # diagD[-NOEFFECTAXIS_INDEX] equals √λ
    if all(centroid == centroid + 0.1 * sigma * diagD[-noEffectAxisIdx] * B[-noEffectAxisIdx]):
        # The coordinate axis std is too low
        restartCriterias["noEffectAxis"] = True

    # Check noEffectCord
    # Stop if adding 0.2-standard deviationos in any single coordinate does not change m
    if any(centroid == centroid + 0.2 * sigma * np.diag(C)):
      # The main axis std has no effect
      restartCriterias["noEffectCoord"] = True
    
    # Check conditionCov
    # Stop if the condition number of the covariance matrix exceeds 10e14
    if cond > 10**14:
      # The condition number is bigger than a threshold
      restartCriterias["conditionCov"] = True


    # Should restar the search
    if any(restartCriterias.values()):
      self.cmaRestart = True
      # Creates a list of restart causes
      restartCauses = [k for k, v in restartCriterias.items() if v]
      print("Restart causes: {}".format(restartCauses))



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
  print("*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*")
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

def printFinalPopulationInfo(status, population, lastFactible, hof, constraintHandling, vtr, maxFe, bestFromLastPopulation):
  print("*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*")
  print("Status: {}".format(status))
  
  print("Value to reach")
  if vtr['individual'] is not None:
    print(vtr['feval'], end=' ')
    vtr['individual'].printIndividual(True, constraintHandling, True)
  else: # Imprime o melhor obtido
    print(maxFe, end=' ')
    hof.printIndividual(True, constraintHandling, True)

  print("Last individual")
  if bestFromLastPopulation is not None:
    print(maxFe, end=' ')
    bestFromLastPopulation.printIndividual(True, constraintHandling, True)
  else:
    print(maxFe, end=' ')
    population.printBest(True, constraintHandling, True)

  print("Hall of fame")
  if hof is not None:
    print(maxFe, end=' ')
    hof.printIndividual(True, constraintHandling, True)
  else:
    print("")

  print("Last factible individual")
  if lastFactible is not None:
    print(maxFe, end=' ')
    lastFactible.printIndividual(True, constraintHandling, True)
  else:
    print("")


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
  else:
    sys.exit("Functiont doesn't have discrete set. Probably should be a continuous one.")
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

def initializeConstraints(function):
  strFunction = str(function)
  truss = lowerBound = upperBound = None
  if strFunction[0] == "1": # Trusses
    truss, lowerBound, upperBound = initializeTruss(function)
    nSize = truss.getDimension()
    g, h, constraintsSize = truss.getNumberConstraints(), 0, truss.getNumberConstraints() + 0
  elif strFunction[0] == "2": # Mechanical Engineering Problems
    if function == 21: # The tension/compression spring design
      nSize = 3
      g, h, constraintsSize = 4, 0, 4
    elif function == 22 or function == 222: # The speed reducer design
      nSize = 7
      g, h, constraintsSize = 11, 0, 11
    elif function == 23: # The welded beam design
      nSize = 4
      g, h, constraintsSize = 5, 0, 5
    elif function == 233: # The welded beam design (Model with 7 constraints and the best objective function value is around 1.7)
      nSize = 4
      g, h, constraintsSize = 7, 0, 7
    elif function == 24: # The pressure vessel design
      nSize = 4
      g, h, constraintsSize = 4, 0, 4
    elif function == 25: # The cantilever beam design
      nSize = 10
      g, h, constraintsSize = 11, 0, 11
    else:
      sys.exit("Function not defined (Initialize Constraints).")
  else:
    sys.exit("Constraints not defined for function.")
  return g, h, constraintsSize, truss, lowerBound, upperBound, nSize # FIXME Verify is none can be returned on trurss, lb and upb

def initiliazeHandleConstraintsParams(function):
  gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize = initializeConstraints(function)
  penaltyCoefficients = [None for i in range(constraintsSize)]
  avgObjFunc = -1
  return gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc

def isInfactible(individual, constraintHandling):
  # Bound constrained prroblems, always factible
  if constraintHandling is None:
    return False
  elif constraintHandling == DEB:
    if individual.violationSum == 0:
      return False
  elif constraintHandling == APM:
    if individual.objectiveFunction[0] == individual.fitness:
      return False
  else:
    sys.exit("Constraint handling method not defined.")
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
    parentsSize = int(4 + 3 * np.log(nSize)) # Population size (λ)
  if centroid is None:
    # centroid = np.array([5.0]*nSize) # objective variables initial points
    centroid = np.random.randn(nSize) # objective variables initial points | gaussian distribution between [0,1]
  if sigma is None:
    # sigma = 5.0 # coordinate wise standard deviation(step-size)
    sigma = 0.5 # coordinate wise standard deviation(step-size)
  if mu is None:
    mu = int(parentsSize / 2) # number of parents selected (selected search points in the population) (µ)
  if rweights is None:
    rweights = "Linear"

  # sys.exit(rweights)
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

  # Strategy parameter setting: Selection (w)
  if rweights == "Superlinear":
    weights = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1))
  elif rweights == "Linear":
    weights = mu + 0.5 - np.arange(1, mu + 1) # muXone recombination weights
  elif rweights == "Equal":
    weights = np.ones(mu)
  else:
    sys.exit("Weights not defined.")
  weights /= np.sum(weights) # normalize recombination weights array
  mueff = 1. / np.sum(weights **2) # vriance-effective size of mu

  # Strategy parameter setting: Adaptation
  cc = 4. / (nSize + 4.) # time constant for cumulation for C (cc)
  cs = (mueff + 2.) / (nSize + mueff + 3.) # t-const for cumulation for sigma control (cσ)
  ccov1 = 2. / ((nSize + 1.3) ** 2 + mueff) # learning rate for rank one update of C
  ccovmu = 2. * (mueff - 2. + 1. / mueff) / ((nSize + 2.) ** 2 + mueff) # and for rank-mu update
  ccovmu = min(1 - ccov1, ccovmu)
  damps = 1. + 2. * max(0, np.sqrt((mueff - 1.) / (nSize + 1.)) - 1.) + cs # # damping for sigma (dσ)

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

# Value to Reach
def valueToReach(function):
  fOptima = None
  if function == 21: # Tension compression spring
    # fOptima = 0.012665 + 1e-8
    # Mezura Montes[07] F_OBJ = 0.01266524482411833 | Project Variables = [11.290483, 0.356692, 0.051688] FOPTIMA
    fOptima = 0.01266524482411833 + 1e-8 # OK
  elif function == 222: # speedReducer
    # fOptima = 2994.471066 + 1e-7
    # Zhang [08] F_OBJ = 2994.471066160767 | Project Variables = [3.5000000000, 0.7000000000, 17, 7.3000000000, 7.7153199115, 3.3502146661, 5.2866544650] # Zhang [08] FOPTIMA
    fOptima = 2994.471066160767 + 1e-8 # OK
  # elif function == 23: # Welded 
  elif function == 233: # Welded Beam Mod(1.7objFunc)
    # fOptima = 1.724852 + 1e-8
    # Wang [10] F_OBJ = 1.7248523110932348 | Project Variables = [0.20572964, 3.47048867, 9.03662391, 0.20572964] FOPTIMA
    fOptima = 1.7248523110932348 + 1e-8 # OK
  elif function == 24: # PressureVesel
    # fOptima = 6059.701660
    fOptima = 6059.701660 + 1e-6 # OK
  else:
    sys.exit("Function not defined (Value To Reach Function).")
  # elif function == 25: # CantileverBeam

  return fOptima
# Generate random indexes that won't repeat on a generation of new offsprings
def populationPick(parentIdx, parentsSize):
  chosenOnes = []
  while len(chosenOnes) < 3:
    idx = np.random.randint(0, parentsSize)
    if idx != parentIdx:
      if idx not in chosenOnes:
        chosenOnes.append(idx)
  return chosenOnes

# Differential evolution
def DE(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case, feval, rweights):
  np.random.seed(seed)
  strFunction = str(function)
  feval = 0
  hof = None
  lastFactible = None
  discreteSet = None
  status = "Initializing"

  if strFunction[0] == "3": # cec2020 bound constrained
    maxFe = defineMaxEval(function, nSize)
    constraintHandling = None

  # Initialize differential evolution parameters
  CR, F, generatedOffspring = deInitParams(function, parentsSize, offspringsSize)

  if case == "discrete":
    discreteSet = getDiscreteCaseList(function)
  # Initialize truss constraints, if necessary
  if constraintHandling:
    gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc = constraintsInitParams(function, constraintHandling)
    # discreteSet = getDiscreteCaseList(function)
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

    # Gets hall of fame individual and last factible one
    hof = parents.hallOfFame(hof, constraintHandling)
    lastFactible = parents.lastFactibleIndividual(lastFactible, constraintHandling)
    
    # Prints best individual of current generation
    parents.printBest(True, constraintHandling, True)
  status = "Finished"
  printFinalPopulationInfo(status, parents, lastFactible, hof, constraintHandling, None)


def CMAES(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case, feval, rweights):
  np.random.seed(seed)
  strFunction = str(function)
  hof = None
  lastFactible = None
  discreteSet = None
  generation = 0
  vtr = {
    'individual': None,
    'feval': None
  }
  status="Initializing"


  
  restartCriterias = {
    'tolFun': False,
    'tolX': False,
    'noEffectAxis': False,
    'noEffectCoord': False,
    'conditionCov': False,
  }
  # if not any(conditions.values()): # Check if should restart

  # stop_causes = [k for k, v in conditions.items() if v]
  # print("Stopped because of condition%s %s" % ((":" if len(stop_causes) == 1 else "s:"), ",".join(stop_causes)))






  # print(AlGORITHM_PARAMS)
  # sys.exit(rweights)

  if strFunction[0] == "3": # cec2020 bound constrained
    maxFe = defineMaxEval(function, nSize)
    constraintHandling = None

  # User defined params (if set to None, uses default)
  # # Search restarted, shouldn't use default population Size
  # if feval != 0: 
  #   parentsSize*=2
  #   centroid = sigma = mu = None
  # else:
  #   parentsSize = centroid = sigma = mu = None
  
  parentsSize = centroid = sigma = mu = None
  # rweights = "equal"
  if case == "discrete":
    discreteSet = getDiscreteCaseList(function)
  if constraintHandling:
    gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize, penaltyCoefficients, avgObjFunc = constraintsInitParams(function, constraintHandling)
  else:
    lowerBound = upperBound = gSize = hSize = truss = None

  
  # Initialize all cmaes params
  parentsSize, mu, centroid, sigma, pc, ps, chiN, C, diagD, B, BD, update_count, weights, mueff, cc, cs, ccov1, ccovmu, damps = cmaInitParams(nSize, parentsSize, centroid, sigma, mu, rweights)

  # tolFunSize = 10 + int(numpy.ceil(30. * nSize / parentsSize)) # last 10 + ⌈30n/λ⌉
  
  # Last 10 + ⌈30n/λ⌉ generations
  populationDeque = deque(maxlen=10 + int(np.ceil(30. * nSize / parentsSize)))

  global AlGORITHM_PARAMS 
  AlGORITHM_PARAMS = {
    "function": function,
    "nSize": nSize,
    "parentsSize": parentsSize,
    "offspringsSize": offspringsSize,
    "seed": seed,
    "maxFe": maxFe,
    "constraintHandling": constraintHandling,
    "case": case,
    "feval": feval,
    "rweights": rweights,
  }

  # Generate initial population
  parents = Population(nSize, parentsSize, function, 1, lowerBound, upperBound, gSize, hSize)

  # valueReached = None

  # CMAESrweights=equal+DEB: T10: 5061.970157299801 | T25: 484.05229151214365 | T60: 311.6048214519769 | T72: 380.0802294286406  
  # CMAES+APM: T10: 5062 | T25: 484 | T60: 313 | T72: 380
  printInitialPopulationInfo("CMA-ES", constraintHandling, function, seed, parentsSize, mu, maxFe, "-", "-", rweights, case)
  if feval > maxFe:
    sys.exit("Maximum number of function evaluations too low.")
  # While fuctions evaluations and restar criterias not met, go on
  while feval < maxFe and not any(restartCriterias.values()):
    # Generation number
    generation += 1
    noEffectAxisIdx = generation % nSize # 
    AlGORITHM_PARAMS["feval"] = feval
    status = "Executing"
    if parents.individuals[0].fitness is not None:
      bestFromLastPopulation = parents.bestIndividual(constraintHandling)
      if not isInfactible(bestFromLastPopulation, constraintHandling) and bestFromLastPopulation.objectiveFunction[0] - valueToReach(function) <= 0:
      # if not isInfactible(bestFromLastPopulation, constraintHandling) and np.round(bestFromLastPopulation.objectiveFunction[0], 6) - valueToReach(function) <= 0:
      # if not isInfactible(bestFromLastPopulation, constraintHandling):
      #   if function === 21:
      #   if bestFromLastPopulation.objectiveFunction[0] <= valueToReach(function):

        if vtr['individual'] is None:
          vtr['individual'] = bestFromLastPopulation
          vtr['feval'] = feval
          # print('Value Reached: {}'.format(feval))
          print('VTR: '.format(vtr))
          # sys.exit("okoko")
          # break
    # Generate new population
    parents.cmaGeneratePopulation(parentsSize, centroid, sigma, BD)

    for individual in parents.individuals:
      # print("Individual no for que buga: {}".format(individual))
      # print("constraintHandling  {}".format(constraintHandling))
      # FIX Code doesnt work when constraintHandling is set to none because np.isnan(None) throws an exception.
      if(np.isnan(individual.objectiveFunction).any() or np.isnan(individual.violationSum).any() or np.isnan(individual.n).any()):
        parents.cmaRestart = True
        printFinalPopulationInfo(status, parents, lastFactible, hof, constraintHandling, bestFromLastPopulation)
        break

    # Check bounds and evaluate
    parents.checkBounds(function, lowerBound, upperBound)
    feval = parents.evaluate(function, feval, truss, case, discreteSet)


    # Handling constraints, if necessary
    if constraintHandling:
      avgObjFunc = parents.handleConstraints(None, constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)
  

    # Sorts population
    parents.sort(None, constraintHandling)

    # Prints generation and its best individual
    print(feval, end=" ")
    parents.printBest(True, constraintHandling, True)

    # Gets bests individuals
    hof = parents.hallOfFame(hof, constraintHandling)
    lastFactible = parents.lastFactibleIndividual(lastFactible, constraintHandling)
    populationDeque.append(parents.bestIndividual(constraintHandling))

    # Update covariance matrix strategy from the population
    try:
      centroid, sigma, pc, ps, B, diagD, C, update_count, BD, cond = parents.cmaUpdateCovarianceMatrix(
        parentsSize, mu, centroid, sigma, weights, mueff, cc, cs, ccov1, ccovmu, damps, pc, 
        ps, B, diagD, C, update_count, chiN, BD
      )
    except:
    # except (RuntimeWarning, RuntimeError, TypeError, np.linalg.LinAlgError):
      break

  
    # Check restart criterias
    # parents.restartCriterias(constraintHandling, restartCriterias, populationDeque, noEffectAxisIdx, centroid, sigma, cond, pc, B, diagD, C)





  # if not any(conditions.values()): # Chegck if should restart


  # Need to restart algorithm?
  if parents.cmaRestart:
    CMAES(AlGORITHM_PARAMS["function"], AlGORITHM_PARAMS["nSize"], AlGORITHM_PARAMS["parentsSize"],
      AlGORITHM_PARAMS["offspringsSize"], AlGORITHM_PARAMS["seed"], AlGORITHM_PARAMS["maxFe"], 
      AlGORITHM_PARAMS["constraintHandling"], AlGORITHM_PARAMS["case"], AlGORITHM_PARAMS["feval"], AlGORITHM_PARAMS["rweights"])

  status = "Finished"
  # Prints final info
  printFinalPopulationInfo(status, parents, lastFactible, hof, constraintHandling, vtr, maxFe, None)
  # print("CPU time used (seconds)")
  # print("666")

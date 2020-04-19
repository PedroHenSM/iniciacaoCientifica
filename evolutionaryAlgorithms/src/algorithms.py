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
import operator as op # For sorting population
from copy import deepcopy

from math import sqrt, log, exp

class Individual(object):
	def __init__(self, n, objectiveFunction, g=None, h=None, violations=None, violationSum=None, fitness=None):
		self.n = n
		self.objectiveFunction = objectiveFunction
		self.g = g
		self.h = h
		self.violations = violations
		self.violationSum = violationSum
		self.fitness = fitness

	def printObjectiveFunction(self):
		for item in self.objectiveFunction:
			print(item, end="	")
		print()

	def printProjectVariables(self):
		for item in self.n:
			print(item, end="	")
		print("")

	def printIndividual(self, boolObjFunc, boolN):
		if boolObjFunc:
			print("{:}\t".format(self.objectiveFunction[0]), end=" ")
		if boolN:
			for n in self.n:
				print("{}\t".format(n), end=" ")
		print()

	# Makes print(individual) a string, not a reference (memory adress)
	def __repr__(self):
		return str(self.__dict__) + "\n"

class Population(object):
	def __init__(self, nSize, popSize, function, objFunctionSize, lowerBound=None, upperBound=None, gSize=None, hSize=None):
		strFunction = str(function)
		self.individuals = []
		for _ in range(popSize):
			n = []
			# objFunc = [None for i in range(objFunctionSize)]
			objFunc = [99999 for i in range(objFunctionSize)]
			g = [99999 for i in range(gSize)]
			h = [99999 for i in range(hSize)]
			violations = [99999 for i in range(gSize + hSize)]
			violationSum = 99999
			for _ in range(nSize):
				if strFunction[0] == "1": # Truss problems
					n.append(np.random.uniform(lowerBound, upperBound))
				elif strFunction[0] == "3": # cec 2020 bound constrained
					n.append(np.random.uniform(-100, 100))
				else:
					sys.exit("Function not defined.")
			self.individuals.append(Individual(n, objFunc, g, h, violations, violationSum))

	def printProjectVariables(self):
		for individual in self.individuals:
			print(individual.n, end="\n")
		print()

	def printObjectiveFunction(self):
		for individual in self.individuals:
			print(individual.objectiveFunction, end="")
		print()

	def evaluate(self, function, fe, truss=None):
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
				# def build_array(a, l, startIdx, size, strFunction):
				# def build_list(l, a, startIdx, size, strFunction):
				# Transfers values from a python list to a C++ array
				populateArray(xArray, individual.n, 0, truss.getDimension())
				valuesList = individual.objectiveFunction + individual.g
				populateArray(valuesArray, valuesList, 0, valuesArraySize)
				# Evaluate population
				truss.evaluation(xArray, valuesArray)
				# Transfers values from a C++ array to a python list
				populateList(individual.n, xArray, 0, truss.getDimension())
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
				populateArray(xArray, individual.n, 0, nSize)
				populateArray(objFuncArray, individual.objectiveFunction, 0, objFuncSize)
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

# TODO Verify way to optime this function (basically the same if its being used three times)
	def selectDE(self, offsprings, generatedOffspring, hasConstraints, constraintHandling=None):
		if hasConstraints:
			if constraintHandling == 1: # Deb
				idx = 0
				for i in range(len(self.individuals)): # Parents size
					bestIdx = idx
					# walks through every n offsprings of each parent
					while idx < generatedOffspring * (i + 1):
						# Picks the best offspring
						if offsprings.individuals[idx].violationSum < offsprings.individuals[bestIdx].violationSum:
							bestIdx = idx
						elif offsprings.individuals[idx].violationSum == offsprings.individuals[bestIdx].violationSum:
							if offsprings.individuals[idx].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
								bestIdx = idx
						idx+=1
					# Best offsprings better than parent, he becomes the parent
					if offsprings.individuals[bestIdx].violationSum < self.individuals[i].violationSum:
						self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
					elif offsprings.individuals[bestIdx].violationSum == self.individuals[i].violationSum:
						if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
							self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
			elif constraintHandling == 2: # APM
				idx = 0
				for i in range(len(self.individuals)): # Parents size
					bestIdx = idx
					# walks through every n offsprings of each parent
					while idx < generatedOffspring * (i + 1):
						# Picks the best offspring
						if offsprings.individuals[idx].fitness < offsprings.individuals[bestIdx].fitness:
							bestIdx = idx
						elif offsprings.individuals[idx].fitness == offsprings.individuals[bestIdx].fitness:
							if offsprings.individuals[idx].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
								bestIdx = idx
						idx+=1
					# Best offsprings better than parent, he becomes the parent
					if offsprings.individuals[bestIdx].fitness < self.individuals[i].fitness:
						self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
					elif offsprings.individuals[bestIdx].fitness == self.individuals[i].fitness:
						if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
							self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])
			else:
				sys.exit("Constraint handling not defined.")
		else:
			idx = 0
			for i in range(len(self.individuals)): # Parents size
				bestIdx = idx
				# walks through every n offsprings of each parent
				while idx < generatedOffspring * (i + 1):
					# Picks the best offspring
					if offsprings.individuals[idx].objectiveFunction[0] < offsprings.individuals[bestIdx].objectiveFunction[0]:
						bestIdx = idx
					idx+=1
				# Best offsprings better than parent, he becomes the parent
				if offsprings.individuals[bestIdx].objectiveFunction[0] < self.individuals[i].objectiveFunction[0]:
					self.individuals[i] = deepcopy(offsprings.individuals[bestIdx])

	def selectCMAES(self, offsprings, hasConstraints):
		if hasConstraints:
			sys.exit("Not implemented.")
		else:
			parentsSize = len(self.individuals)
			offspringsSize = len(offsprings.individuals)
			popSize = parentsSize + offspringsSize
			population = Population(len(self.individuals[0].n), popSize, 31, 1)
			for i in range(parentsSize + offspringsSize):
				if i < parentsSize:
					population.individuals[i] = deepcopy(self.individuals[i])
				else:
					population.individuals[i] = deepcopy(offsprings.individuals[i-parentsSize])

			population.sort(None)
			for i in range(parentsSize):
				self.individuals[i] = deepcopy(population.individuals[i])

	def checkBounds(self, function, lowerBound=None, upperBound=None):
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

	def sort(self, offsprings=None, constraintHandling=None):
		if constraintHandling is None: # Bound constraint problem
			self.individuals.sort(key=op.attrgetter("objectiveFunction"))
			if offsprings is not None:
				offsprings.individuals.sort(key=op.attrgetter("objectiveFunction"))
		elif constraintHandling == 1: # Deb
			self.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
			if offsprings is not None:
				offsprings.individuals.sort(key=op.attrgetter("violationSum", "objectiveFunction"))
		elif constraintHandling == 2: # APM
			self.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))
			if offsprings is not None:
				offsprings.individuals.sort(key=op.attrgetter("fitness", "objectiveFunction"))
		else:
			sys.exit("Constraint handling not defined.")

	def bestIndividual(self, constraintHandling=None):
		# Assume the first individual is the best one
		best = deepcopy(self.individuals[0])
		for i in range(1, len(self.individuals)):
			if constraintHandling == 1: # Deb
				if self.individuals[i].violationSum < best.violationSum:
					best = deepcopy(self.individuals[i])
				elif self.individuals[i].violationSum == best.violationSum:
					if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
						best = deepcopy(self.individuals[i])
			elif constraintHandling == 2: # APM
				if self.individuals[i].fitness < best.fitness:
					best = deepcopy(self.individuals[i])
				elif self.individuals[i].fitness == best.fitness:
					if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
						best = deepcopy(self.individuals[i])
			else:
				if self.individuals[i].objectiveFunction[0] < best.objectiveFunction[0]:
					best = deepcopy(self.individuals[i])
		return best

	def hallOfFame(self, hof):
		currentBest = self.bestIndividual()
		if hof is None or currentBest.objectiveFunction[0] < hof.objectiveFunction[0]:
			hof = deepcopy(currentBest)
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
					if constraintHandling == 1: # Deb
						individual.violations[i] = individual.h[idxH]
					elif constraintHandling == 2: # APM
						individual.violations[i] = np.abs(individual.h[idxH]) - 0.0001
					idxH+=1
			if constraintHandling == 1: # Deb
				# Only sums positives values. Negative violations aren't summed
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
			denominator = denominator + sumViolation[i] *  sumViolation[i]
		
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

	def handleConstraints(self, constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc):
		self.uniteConstraints(constraintHandling)
		avgObjFunc = self.calculatePenaltyCoefficients(constraintsSize, penaltyCoefficients, avgObjFunc)
		self.calculateAllFitness(constraintsSize, penaltyCoefficients, avgObjFunc)
		return avgObjFunc

	def printBest(self, boolObjFunc, boolN, constraintHandling=None):
		best = self.bestIndividual(constraintHandling)
		if boolObjFunc:
			print("{:}\t".format(best.objectiveFunction[0]), end=" ")
		if boolN:
			for n in best.n:
				print("{}\t".format(n), end=" ")
		print()

	def moveArzToPop(self, list2d):
		for i in range(len(self.individuals)): # Lambda_
			self.individuals[i].n = deepcopy(list2d[i])

	def selectMuIndividuals(self, mu):
		population = []
		for i in range(mu):
			population.append(self.individuals[i].n)
		return population

	def cmaPopulationGenerate(self, parentsSize, centroid, sigma, BD):
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
				+ sqrt(cs * (2 - cs) * mueff) / sigma \
				* np.dot(B, (1. / diagD) *
										np.dot(B.T, c_diff))

		hsig = float((np.linalg.norm(ps) /
									sqrt(1. - (1. - cs) ** (2. * (update_count + 1.))) / chiN <
									(1.4 + 2. / (nSize + 1.))))

		update_count += 1

		pc = (1 - cc) * pc + hsig \
				* sqrt(cc * (2 - cc) * mueff) / sigma \
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

	# Makes print(population) a string, not a reference (memory adress)
	def __repr__(self):
		return str(self.__dict__) + "\n"

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

def cmaInitParams(nSize, parentsSize, centroid, sigma, mu, rweights):
	# User defined params
	if parentsSize is None:
		parentsSize = int(4 + 3 * log(nSize)) # Population size
	if centroid is None:
		centroid = np.array([5.0]*nSize) # objective variables initial points
	if sigma is None:
		sigma = 5.0 # coordinate wise standard deviation(step-size)
	if mu is None:
		mu = int(parentsSize / 2) # number of parents selected (selected search points in the population)
	if rweights is None:
		rweights = "linear"

	# Initiliaze dynamic (internal) strategy parameters and constants
	pc = np.zeros(nSize) # evolution paths for C
	ps = np.zeros(nSize) # evolutions paths for sigma
	chiN = sqrt(nSize) * (1 - 1. / (4. * nSize) + 1. / (21. * nSize ** 2)) # expectation of ||N(0,I)|| == norm(randn(N,1))
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
		weights = log(mu + 0.5) - np.log(np.arange(1, mu + 1))
	elif rweights == "linear":
		weights = log(mu + 0.5) - np.log(np.arange(1, mu + 1)) # muXone recombination weights
	elif rweights == "equal":
		weights = np.ones(mu)
	weights /= sum(weights) # normalize recombination weights array
	mueff = 1. / sum(weights **2) # vriance-effective size of mu

	# Strategy parameter setting: Adaptation
	cc = 4. / (nSize + 4.) # time constant for cumulation for C
	cs = (mueff + 2.) / (nSize + mueff + 3.) # t-const for cumulation for sigma control
	ccov1 = 2. / ((nSize + 1.3) ** 2 + mueff) # learning rate for rank one update of C
	ccovmu = 2. * (mueff - 2. + 1. / mueff) / ((nSize + 2.) ** 2 + mueff) # and for rank-mu update
	ccovmu = min(1 - ccov1, ccovmu)
	damps = 1. + 2. * max(0, sqrt((mueff - 1.) / (nSize + 1.)) - 1.) + cs # # damping for sigma

	return parentsSize, mu, centroid, sigma, pc, ps, chiN, C, diagD, B, BD, update_count, weights, mueff, cc, cs, ccov1, ccovmu, damps

# Python function to pass values of a python list to a C++ (swig) array
def populateArray(a, l, startIdx, size):
	for i in range(size - startIdx):
		# Sets the arr[startIdx] with value from list[i]
		utils.doubleArray_setitem(a, startIdx, l[i])
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
			sys.exit("Dimension not defined for cec2020 bounded constrained functions")
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

def initializeConstraints(function, nSize):
	strFunction = str(function)
	truss = lowerBound = upperBound = None
	if strFunction[0] == "1": # Trusses
		truss, lowerBound, upperBound = initializeTruss(function)
		nSize =  truss.getDimension()
		g, h, constraintsSize = truss.getNumberConstraints(), 0, truss.getNumberConstraints() + 0
	else:
		sys.exit("Constraints not defined for function.")
	return g, h, constraintsSize, truss, lowerBound, upperBound, nSize

# Differential evolution
def DE(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling=None):
	np.random.seed(seed)
	strFunction = str(function)
	hasConstraints = False
	if strFunction[0] == "3": # cec2020 bound constrained
		maxFe = defineMaxEval(function, nSize)

	# Define DE parameters
	CR = 0.9
	F = 0.6
	feval = 0
	hof = None
	generatedOffspring = int(offspringsSize / parentsSize)
	if strFunction[0] == "1":
		hasConstraints = True
	truss = lowerBound = upperBound = None
	if hasConstraints:
		gSize, hSize, constraintsSize, truss, lowerBound, upperBound, nSize = initializeConstraints(function, nSize)
		penaltyCoefficients = [None for i in range(constraintsSize)]
		avgObjFunc = -1
		parents = Population(nSize, parentsSize, function, 1, lowerBound, upperBound, gSize, hSize)
		offsprings = Population(nSize, offspringsSize, function, 1, lowerBound, upperBound, gSize, hSize)
		feval = parents.evaluate(function, feval, truss)
		feval = offsprings.evaluate(function, feval, truss)
		# Starts constraints handling
		if constraintHandling == 1: # Deb penalization
			avgObjFunc = offsprings.handleConstraints(constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)
			offsprings.sort(constraintHandling=constraintHandling)
			parents.selectDE(offsprings, generatedOffspring, True, 1) # TODO Remove "hasConstraints" parameter?
		elif constraintHandling == 2: # APM
			avgObjFunc = parents.handleConstraints(constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)
		else:
			sys.exit("Constraint handling method not defined.")
	else:
		parents = Population(nSize, parentsSize, function, 1)
		offsprings = Population(nSize, offspringsSize, function, 1)
		feval = parents.evaluate(function, feval, truss)

	while feval < maxFe:
		offspringIdx = 0
		for i in range(parentsSize):
			for _ in range(generatedOffspring):
				chosenOnesIdxs = populationPick(i, parentsSize)
				R = np.random.randint(0, nSize)  # Random index
				for j in range(nSize):
					Ri = np.random.rand()  # Generates random number between (0,1)
					if Ri < CR or j == R:
						offsprings.individuals[offspringIdx].n[j] = parents.individuals[chosenOnesIdxs[0]].n[j] + F * (parents.individuals[chosenOnesIdxs[1]].n[j] - parents.individuals[chosenOnesIdxs[2]].n[j])
					else:
						offsprings.individuals[offspringIdx].n[j] = parents.individuals[i].n[j]
				offspringIdx = offspringIdx + 1

		# Check bounds and evaluate offsprings
		offsprings.checkBounds(function, lowerBound, upperBound)
		feval = offsprings.evaluate(function, feval, truss)
		
		if hasConstraints:
			if constraintHandling == 1: # Deb
				avgObjFunc = offsprings.handleConstraints(constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)
			elif constraintHandling == 2: # APM
				avgObjFunc = offsprings.handleConstraints(constraintHandling, constraintsSize, penaltyCoefficients, avgObjFunc)

		# Selects best offsprings and put them in parents
		parents.selectDE(offsprings, generatedOffspring, hasConstraints, constraintHandling)

		# Gets hall of fame individual
		hof = parents.hallOfFame(hof)
		# Prints the best individual of the generation
	# parents.printBest(True, False, constraintHandling)
	hof.printIndividual(True, False)


def CMAESOld(function, nSize, parentsSize, offspringsSize, seed, maxFe):
	np.random.seed(seed)
	strFunction = str(function)
	hasConstraints = False
	if strFunction[0] == "3": # cec2020 bound constrained
		maxFe = defineMaxEval(function, nSize)

	feval = 0
	generatedOffspring = int(offspringsSize / parentsSize)

	# User defined parameters
	sigma = 0.5 # coordinate wise standard deviation (step-size)
	xmean = np.random.randn(nSize)  # objective variables initial point

	# λ ≥ 2, population size, sample size, number of offspring, see (5).
	# µ ≤ λ parent number, number of (positively) selected search points in the population, number of strictly positive recombination weights, see (6).
	
	# λ = tamanho da populacao
	# µ = filhos escolhidos

	# # Strategy parameters setting: Selection
	parentsSize = 4 + np.floor(3 * np.log(nSize))  # population size (λ), offsprings number (same as parents)
	parentsSize = int(parentsSize)
	# parentsSize = nSize * 4 # population size, offsprings number (same as parents)
	offspringsSize = 4 * parentsSize 



	if hasConstraints:
		sys.exit("Not implemented.")
	else:
		parents = Population(nSize, parentsSize, function, 1)
		offsprings = Population(nSize, offspringsSize, function, 1)
		feval = parents.evaluate(function, feval)

	# Strategy parameters setting: Selection

	# mu = µ = offspringsSize

	mu = parentsSize / 2   # number of parents selected (selected search points in the population)
	muList = [i + 1 for i in range(int(mu))]
	weights = np.log(mu+1/2)-np.log(muList).conj().T  # muXone recombination weights
	mu = np.floor(mu) # number of parents selected (µ)(selected search points in the population)
	weights = weights/np.sum(weights) # normalize recombination weights array
	mueff = np.power(np.sum(weights), 2) / np.sum(np.power(weights, 2)) # variance-effective size of mu

	# Strategy parameter setting: Adaptation
	cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
	cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
	c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
	cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
	damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

	# Initiliaze dynamic (internal) strategy parameters and constants
	pc = np.zeros(nSize) # evolutions paths for C
	ps = np.zeros(nSize)  # evolutions paths for sigma
	B = np.eye(nSize)  # B defines de coordinate system
	D = np.eye(nSize)  # diagonal matrix D defines the scaling
	C = B @ D @ (B@D).conj().T # covariance matrix
	eigenval = 0  # B and D update at counteval == 0
	chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))


	hof = None

	deg = 0
	while feval < maxFe:
		arzAuxList = []
		arxAuxList = []
		for i in range(offspringsSize):
			arz = np.random.randn(nSize)  # standard normally distributed vector  (38)
			# print("Arz: {}".format(arz))
			arx = xmean + sigma * (B @ D @ arz)  # add mutation N(m,sigma²C) (40)
			arzAuxList.append(arz)
			arxAuxList.append(arx)
			for j in range(nSize):  # Copies arx to individual[].n TODO This can be done as offsprings.individuals[k].n = [i for i in arx] ?!
				offsprings.individuals[i].arz[j] = arz[j]
				offsprings.individuals[i].n[j] = arx[j]
				complexx = np.iscomplexobj(offsprings.individuals[i].n)
		
		if complexx:
			hof.printIndividual(True, False)
			sys.exit("Complex number in individual.")
		
		arz = np.stack(arzAuxList)
		arx = np.stack(arxAuxList)

		# Check bounds and evaluate offsprings
		offsprings.checkBounds(function)
		feval = offsprings.evaluate(function, feval)

		# parents.sort(offsprings)
		if hasConstraints:
			sys.exit("Not implemented.")
		else:
			# Picks best individuals among parents and offsprings and put it on parents
			parents.selectCMAES(offsprings, False)
		
		arz = arz.T
		arx = arx.T

		for j in range(parentsSize):
			for i in range(nSize):
				arx[i][j] = parents.individuals[j].n[i]
				arz[i][j] = parents.individuals[j].arz[i] 

		# Individuals already stored in columns, no need to tranpose
		muBestX = np.delete(arx, np.s_[int(mu):], 1)  # Remove columns after mu index (Select mu individuals)
		muBestZ = np.delete(arz, np.s_[int(mu):], 1)  # Remove columns after mu index (Select mu individuals)
		xmean = np.matmul(muBestX, weights)  # xmean is one array with nSize positions. Recombination Eq. 42
		zmean = np.matmul(muBestZ, weights)  # zmeanis one array with nSize positions. == Dˆ-1*B’*(xmean-xold)/sigma

		# Cumulation: Updatte evolution paths
		ps = (1-cs)*ps + (np.sqrt(cs*(2-cs)*mueff)) * B@zmean  # Eq. 43
		hsig = True if np.linalg.norm(ps) / np.sqrt(1-np.power((1-cs), (2*feval/parentsSize)))/chinN < 1.4 + 2/(nSize + 1) else False
		pc = (1-cc)*pc + hsig * np.sqrt(cc*(2-cc)*mueff) * B@D@zmean # Eq. 45
		
		# Adapt covariance matrix C. Eq. 47
		# C = (1-c1-cmu) * C + c1 * (np.outer(pc, pc) + (1-hsig) * cc*(2-cc) * C) + cmu * (B@D@muBestZ) @ np.diag(weights) @ (B@D@muBestZ).conj().T # Eq. 47
		C = ((1-c1-cmu) * C + # regard old matrix
			c1 * (np.outer(pc, pc) + # plus rank on update
			(1-hsig) * cc*(2-cc) * C) + # minor correction
			cmu *  # plus rank mu update
			(B@D@muBestZ) @ np.diag(weights) @ (B@D@muBestZ).conj().T)


		#  Adapt step-size sigma
		sigma = sigma * np.exp((cs/damps) * (np.linalg.norm(ps)/chinN - 1)) # Eq. 44

		# Gets the hall of fame individual
		if (hof is None):
			hof = parents.bestIndividual()
		else:
			bestFromActualGen = parents.bestIndividual()
			if bestFromActualGen.objectiveFunction[0] < hof.objectiveFunction[0]:
				hof = deepcopy(bestFromActualGen)

		if feval - eigenval > parentsSize / (c1 + cmu) / nSize / 10:  # to achieve 0(N^2)
			eigenval = feval

			C = np.triu(C) + np.triu(C, 1).conj().T  # enforce symmetry
			D, B = np.linalg.eig(C)  # eigen decomposition, B == normalized eigenvector 

		if(np.any(D<=0)): # degenerates
			deg = deg + 1
			print("Degenerou {} vezes".format(deg))

			# User defined parameters
			sigma = 0.5 # coordinate wise standard deviation (step-size)
			xmean = np.random.randn(nSize)  # objective variables initial point

			# Strategy parameter setting: Adaptation
			cc = (4+mueff / nSize) / (nSize+4 + 2*mueff/nSize) # time constant for cumulation for C
			cs = (mueff+2) / (nSize+mueff+5)  # t-const for cumulation for sigma control
			c1 = 2 / (np.power(nSize + 1.3, 2) + mueff)  # learning rate for rank one update of C
			cmu = np.minimum(1 - c1, 2 * (mueff - 2 + 1/mueff)/(np.power(nSize+2, 2) + 2*mueff/2))  # and for rank-mu update
			damps = 1 + 2*np.maximum(0, np.sqrt((mueff -1) / (nSize + 1)) -1 ) + cs  # damping for sigma

			# Initiliaze dynamic (internal) strategy parameters and constants
			pc = np.zeros(nSize) # evolutions paths for C
			ps = np.zeros(nSize)  # evolutions paths for sigma
			B = np.eye(nSize)  # B defines de coordinate system
			D = np.eye(nSize)  # diagonal matrix D defines the scaling
			C = B @ D @ (B@D).conj().T # covariance matrix
			# eigenval = 0  # B and D update at counteval == 0 # TODO Verify is eigenval need to be 0 again.
			chinN = nSize**(0.5) * (1-1/(4*nSize)+1 / (21*np.power(nSize, 2)))  # expectation of ||N(0,I)|| == norm(randn(N,1))
		else:
			D = np.diag(np.sqrt(D)) # D containts standard deviations now
	hof.printIndividual(True, False)

def CMAES(function, nSize, parentsSize, offspringsSize, seed, maxFe):
	np.random.seed(seed)
	strFunction = str(function)
	hasConstraints = False
	if strFunction[0] == "3": # cec2020 bound constrained
		maxFe = defineMaxEval(function, nSize)

	feval = 0
	# User defined params (if set to None, uses default)
	parentsSize = centroid = sigma = mu = rweights = None
	rweights = "equal"
	# Init all cmaes params
	parentsSize, mu, centroid, sigma, pc, ps, chiN, C, diagD, B, BD, update_count, weights, mueff, cc, cs, ccov1, ccovmu, damps = cmaInitParams(nSize, parentsSize, centroid, sigma, mu, rweights)

	hof = None
	if hasConstraints:
		sys.exit("Not implemented.")
	else:
		parents = Population(nSize, parentsSize, function, 1) # TODO Melhorar isso (Gerando 2x)

	while feval < maxFe:
		# print(feval)
		# Generate new population
		parents.cmaPopulationGenerate(parentsSize, centroid, sigma, BD)
		# Evaluate and sort
		feval = parents.evaluate(function, feval)
		parents.sort()
		# Gets the 'hall of fame" individual
		hof = parents.hallOfFame(hof)
		# NOTE Update covariance matrix strategy from the population
		centroid, sigma, pc, ps, B, diagD, C, update_count, BD = parents.cmaUpdateCovarianceMatrix(
			parentsSize, mu, centroid, sigma, weights, mueff, cc, cs, ccov1, ccovmu, damps, pc, 
			ps, B, diagD, C, update_count, chiN, BD
		)
		# parents.printBest(True, True)
	hof.printIndividual(True, True)
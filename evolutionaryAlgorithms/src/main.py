import sys
import algorithms
import argparse


import numpy as np
from math import sqrt, log, exp

def execAlgorithm(algorithm, function, nSize, parentsSize, offspringsSize, seed, maxFe):
	if algorithm == "GA":
		sys.exit("Not implemented.")
	elif algorithm == "DE":
		algorithms.DE(function, nSize, parentsSize, offspringsSize, seed, maxFe)
	elif algorithm == "ES":
		sys.exit("Not implemented.")
	elif algorithm == "CMAES":
		algorithms.CMAES(function, nSize, parentsSize, offspringsSize, seed, maxFe)
	else:
		sys.exit("Algorithm not found.")


def menu():
	parser = argparse.ArgumentParser(description="Evolutionary Algorithms")
	parser.add_argument("--algorithm", "-a", type=str, default="CMAES", help="Algorithm to be used (GA, ES, DE or CMAES)")
	parser.add_argument("--function", "-f", type=int, default=31, help="Truss to be solved (10, 25, 60, 72 or 942 bars). "
						"For the truss problem, the first digit must be 2, followed by the number of the bars in the problem. "
						"Example: 225, is for the truss of 25 bars")
	parser.add_argument("--seed", "-s", type=int, default=1, help="Seed to be used")
	# parser.add_argument("--penaltyMethod", "-p", type=int, default=1, help="Penalty method to be used. 1 for Deb Penalty or 2 for APM")
	parser.add_argument("--parentsSize", "-u", type=int, default=50, help="µ is the parental population size")  # u from µ (mi) | µ ≈ λ/4
	parser.add_argument("--nSize", "-n", type=int, default=5, help="Search space dimension")
	parser.add_argument("--offspringsSize", "-l", type=int, default=50, help="λ is number of offsprings, offsprings population size")  # l from λ (lambda) | µ ≈ λ/4
	parser.add_argument("--maxFe", "-m", type=int, default=10000, help="The max number of functions evaluations")
	# parser.add_argument("--crossoverProb", "-c", type=int, default=100, help="The crossover probability [0,100]")
	# parser.add_argument("--esType", "-e", type=int, default=0, help="The type of ES. 0 for ES(µ + λ) or 1 for ES(µ , λ)")
	# parser.add_argument("--globalSigma", "-g", type=int, default=0, help="If the σ parameter is global or not. 1 for global σ or 0 if not")
	# parser.add_argument("--windowSize", "-w", type=int, default=5, help="Size of the window for updating gaussian model")
	args = parser.parse_args()
	# CEC20 Bound Constrained
	# F1-F5 & F8-F10 : D = 5, 10, 15, 20
	# F6 & F7 : D = 10, 15,
	execAlgorithm(args.algorithm, args.function, args.nSize, args.parentsSize, args.offspringsSize, args.seed, args.maxFe)
	# algorithm(args.algorithm, args.function, args.seed, args.penaltyMethod, args.parentsSize, args.nSize, args.offspringsSize, args.maxFE, args.crossoverProb, args.esType, args.globalSigma, args.windowSize)
	print(args)



	
if __name__ == '__main__':
	# l = [10, 20, 30, 40, 50]
	# for idx, item in enumerate(l):
	# 	print("Item {}: {}".format(idx, item))
	# for idx, item in enumerate(l):
	# 		l[idx] = -9
	# print(l)
	# p1 = algorithms.Population(5, 3, 31, 1)
	# p2 = algorithms.Population(5, 5, 31, 1)
	# p1.individuals[0].objectiveFunction[0] = 10
	# p1.individuals[1].objectiveFunction[0] = 20
	# p1.individuals[2].objectiveFunction[0] = 30

	# p2.individuals[0].objectiveFunction[0] = 15
	# p2.individuals[1].objectiveFunction[0] = -20
	# p2.individuals[2].objectiveFunction[0] = 34
	# p2.individuals[3].objectiveFunction[0] = 50
	# p2.individuals[4].objectiveFunction[0] = -4
	# print(p1)
	# print(p2)
	# p1.selectCMAES(p2,False)
	# print(p1)




	menu()







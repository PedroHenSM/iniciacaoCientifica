import sys
import algorithms
import argparse


import numpy as np
import timeit
from math import sqrt, log, exp

def execAlgorithm(algorithm, function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case):
  if algorithm == "GA":
    sys.exit("Not implemented.")
  elif algorithm == "DE":
    startTime = timeit.default_timer()
    algorithms.DE(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case, 0)
    endTime = timeit.default_timer()
    elapsedTime = endTime - startTime
    print("CPU time used (seconds)\n{}".format(elapsedTime))

    # algorithms.DE(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling)
    # algorithms.DE(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling)
  elif algorithm == "ES":
    sys.exit("Not implemented.")
  elif algorithm == "CMAES":
    startTime = timeit.default_timer()
    algorithms.CMAES(function, nSize, parentsSize, offspringsSize, seed, maxFe, constraintHandling, case, 0)
    endTime = timeit.default_timer()
    elapsedTime = endTime - startTime
    print("CPU time used (seconds)\n{}".format(elapsedTime))
  else:
    sys.exit("Algorithm not found.")


def menu():
  parser = argparse.ArgumentParser(description="Evolutionary Algorithms")
  parser.add_argument("--algorithm", "-a", type=str, default="CMAES", help="Algorithm to be used (GA, ES, DE or CMAES)")
  parser.add_argument("--function", "-f", type=int, default=31, help="Truss to be solved (10, 25, 60, 72 or 942 bars). "
            "For the truss problem, the first digit must be 2, followed by the number of the bars in the problem. "
            "Example: 225, is for the truss of 25 bars")
  parser.add_argument("--seed", "-s", type=int, default=1, help="Seed to be used")
  parser.add_argument("--constraintHandling", "-p", type=str, default="DEB", help="Constraint handling method to be used (DEB or APM)")
  parser.add_argument("--parentsSize", "-u", type=int, default=50, help="µ is the parental population size")  # u from µ (mi) | µ ≈ λ/4
  parser.add_argument("--nSize", "-n", type=int, default=5, help="Search space dimension")
  parser.add_argument("--offspringsSize", "-l", type=int, default=50, help="λ is number of offsprings, offsprings population size")  # l from λ (lambda) | µ ≈ λ/4
  parser.add_argument("--maxFe", "-m", type=int, default=15000, help="The max number of functions evaluations")
  parser.add_argument("--case", "-c", type=str, default="continuous", help="Discrete or continuous design variables")
  # parser.add_argument("--crossoverProb", "-c", type=int, default=100, help="The crossover probability [0,100]")
  # parser.add_argument("--esType", "-e", type=int, default=0, help="The type of ES. 0 for ES(µ + λ) or 1 for ES(µ , λ)")
  # parser.add_argument("--globalSigma", "-g", type=int, default=0, help="If the σ parameter is global or not. 1 for global σ or 0 if not")
  # parser.add_argument("--windowSize", "-w", type=int, default=5, help="Size of the window for updating gaussian model")
  args = parser.parse_args()
  # CEC20 Bound Constrained
  # F1-F5 & F8-F10 : D = 5, 10, 15, 20
  # F6 & F7 : D = 10, 15,
  execAlgorithm(args.algorithm, args.function, args.nSize, args.parentsSize, args.offspringsSize, args.seed, args.maxFe, args.constraintHandling, args.case)
  # print(args)



  
if __name__ == '__main__':
  # 15 000 maxFe
  # 30 seeds
  # t10 and t72 used 10 000 maxFe
  # t942 used 150 000 maxFe
  # SMDE with KNN is the best one - page 24 Krempser



  menu()







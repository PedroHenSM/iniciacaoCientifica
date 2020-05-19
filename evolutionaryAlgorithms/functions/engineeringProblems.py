import numpy as np
import sys

def f01(x, objFunc, g, h): # The tension/compression spring design
  objectiveFunction = (x[0] + 2) * x[1] * np.power(x[2], 2) # Volume
  g[0] = 1 - ((np.power(x[1], 3)*x[0])/(71785*np.power(x[2],4))) # <= 0
  g[1] = ((4*np.power(x[1],2) - x[2]*x[1])/(12566(x[1] * np.power(x[2],3) - np.power(x[2],4)))) + (1/(5108 * np.power(x[2], 2))) - 1 # <= 0
  g[2] = 1 - ((140.45) * x[2] / (np.power(x[1], 2) * x[0])) # <= 0
  g[3] = ((x[1] + x[2]) / 1.5) - 1 # <= 0

  return objectiveFunction

def f02(x, objFunc, g, h): # The speed reducer design
  objectiveFunction = 0.7854 * x[0] * np.power(x[1], 2) * (3.3333 * np.power(x[2], 2) + 14.9334 * x[2] - 43.0934) - 1.508 * x[0] * (np.power(x[5], 2) + np.power(x[6], 2)) + 7.4777 * (np.power(x[5], 3) + np.power(x[6], 3)) + 0.7854 * (x[3] * np.power(x[5], 2) + x[4] * np.power(x[6],2)) # Weight
  g[0] = 27 * np.float_power(x[0], -1) * np.float_power(x[1], -2) * np.float_power(x[2], -1) -1 # <= 1
  g[1] = 397.5 * np.float_power(x[0], -1) * np.float_power(x[1], -2) * np.float_power(x[2], -2) -1 # <= 1
  g[2] = 1.93 * np.float_power(x[1], -1) * np.float_power(x[2], -1) * np.power(x[3], 3) * np.float_power(x[5], -4) -1 # <= 1
  g[3] = 1.93 * np.float_power(x[1], -1) * np.float_power(x[2], -1) * np.power(x[4], 3) * np.float_power(x[6], -4) -1 # <= 1
  g[4] = (1/(0.1 * np.power(x[5], 3))) * np.power((np.power(((745 * x[3]) / (x[1] * x[2])),2) + 16.9*10e6 ), 0.5) -1100 # <=1100
  g[5] = (1/(0.1 * np.power(x[6], 3))) * np.power((np.power(((745 * x[4]) / (x[1] * x[2])),2) + 157.5*10e6 ), 0.5) -850 # <=850
  g[6] = x[1] * x[2] - 40 # <= 40
  g[7] = 5 - (x[0] / x[1]) # >= 5
  g[8] = x[0] / x[1] -12 # <= 12
  g[9] = (1.5 * x[5] + 1.9) * np.float_power(x[3], -1) -1 # <= 1
  g[10] = (1.1 * x[6] + 1.9) * np.float_power(x[4], -1) -1 # <= 1
  
  return objectiveFunction

def f03(x, objFunc, g, h): # The welded beam design
  # h -> x[0] | l -> x[1] | t -> x[2] | b -> x[3]
  objectiveFunction = 1.10471 * np.power(x[0], 2) * x[1] + 0.04811 * x[2] * x[3](14.0 + x[1]) # Cost
  g[0] = 13600 - np.sqrt(np.power(tau1, 2) + np.power(tau2, 2) + x[1] * tau1 * tau2 / alpha) # >= 0
  g[1] = 30000 - 504000/(np.power(x[2], 2) * x[3]) # >= 0
  g[2] = x[3] - x[0] # >= 0
  g[3] = pc - 6000 # >= 0
  g[4] = 0.25 - 2.1952/(np.power(x[2],3) * x[3]) # >= 0

  alpha = np.sqrt(0.25 * (np.power(x[1], 2) + np.power((x[0] + x[2]), 2)))
  tau1 = (6000) / (np.sqrt(2) * x[0] * x[1])
  tau2 = ((6000 * (14 + 0.5 * x[1]) * alpha) / (2 * (0.707 * x[0] * x[1] * (np.power(x[1], 2) / 12 + 0.25 * np.power((x[0] + x[2]),2)))))
  pc = 64746.022(1 - 0.0282346 * x[2]) * x[2] * np.power(x[3], 3)

  return objectiveFunction

def f04(x, objFunc, g, h): # The pressure vessel design
  # Ts -> x[0] | Th -> x[1] | R -> x[2] | L -> x[3]
  w = 0, 6224 * x[0] * x[1] * x[2] + 1.7781 * x[1] * np.power(x[2], 2) + 3.1661 * np.power(x[0], 2) * x[3] + 19.84 * np.power(x[0], 2) * x[2] # Weight
  g[0] = x[0] - 0.0193 * x[2] # >= 0
  g[1] = x[1] - 0.00954 * x[2] # >= 0
  g[2] = np.pi * np.power(x[2], 2) * x[3] + 4/(3 * np.pi * np.power(x[2], 3)) - 1296000 # >= 0
  g[3] = -x[3] + 240 # >= 0

  return objectiveFunction

def f05(x, objFunc, g, h): # The cantilever beeam design
  # Hi -> x[0] | Bi -> x[1]
  objectiveFunction = 0 # Volume
  for i in range(5):
    v += x[0] * x[1]
  objectiveFunction = objectiveFunction * 100

  g[0] = g[1] = g[2] = g[3] = g[4] = sigmai - 14000 # <= 14000 # FIXME Verify value of sigmai

  g[5] = g[6] = g[7] = g[8] = g[9] = x[0] / x[1] - 20 # <= 20
  g[10] = delta - 2.7 # <= 2.7

  return objectiveFunction

def executeFunction(function, x, objFunc, g, h):
  if function == 21:
    objFunc, g, h = f01(x, objFunc, g, h)
  elif function == 22:
    objFunc, g, h = f02(x, objFunc, g, h)
  elif function == 23:
    objFunc, g, h = f03(x, objFunc, g, h)
  elif function == 24:
    objFunc, g, h = f04(x, objFunc, g, h)
  elif function == 25:
    objFunc, g, h = f05(x, objFunc, g, h)
  else:
    sys.exit("Function not defined.")
  return objFunc, g, h

# if __name__ == '__main__':
#   executeFunction()










# *---*---*---*---*---*---*---*---*---* Pandas Personal Helper ---*---*---*---*---*---*---*---*---*---*---*

# g(x) <= V -> V - g(x) <= 0
# g(x) >= V -> g(x) - V >= 0 

# def rosen():
# 	if function == 91:  # RosenbrockFunction
# 	sumRosen = 0
# 	for k in range(nSize - 1):
# 		sumRosen = sumRosen + 100*np.power((self.individuals[i].n[k+1] - np.power(self.individuals[i].n[k], 2)), 2) + np.power((1 - self.individuals[i].n[k]), 2)
# 	self.individuals[i].objectiveFunction[0] = sumRosen
import numpy as np
import sys

# Helper function
def find_nearest(array, value):
  array = np.asarray(array)
  idx = (np.abs(array - value)).argmin()
  return array[idx]

def tensionCompressionSpring(x, objFunc, g, h): # The tension/compression spring design | 36k maxEval
  objectiveFunction = (x[0] + 2) * x[1] * np.power(x[2], 2) # Volume
  g[0] = 1 - ((np.power(x[1], 3)*x[0])/(71785*np.power(x[2],4))) # <= 0
  g[1] = ((4*np.power(x[1],2) - x[2]*x[1])/(12566 * (x[1] * np.power(x[2],3) - np.power(x[2],4)))) + (1/(5108 * np.power(x[2], 2))) - 1 # <= 0
  g[2] = 1 - ((140.45) * x[2] / (np.power(x[1], 2) * x[0])) # <= 0
  g[3] = ((x[1] + x[2]) / 1.5) - 1 # <= 0

  return objectiveFunction, g, h

def speedReducer(x, objFunc, g, h): # The speed reducer design | 36k maxEval
  # x =[3.5, 0.7, 17, 7.3, 7.8, 3.350215, 5.286683] Best solution found by bernardino (08)
  # x[2] = np.round(x[2]) # This variable is integer

  intx2 = np.round(x[2])
  objectiveFunction = 0.7854 * x[0] * np.power(x[1], 2) * (3.3333 * np.power(intx2, 2) + 14.9334 * intx2 - 43.0934) - 1.508 * x[0] * (np.power(x[5], 2) + np.power(x[6], 2)) + 7.4777 * (np.power(x[5], 3) + np.power(x[6], 3)) + 0.7854 * (x[3] * np.power(x[5], 2) + x[4] * np.power(x[6],2)) # Weight
  g[0] = 27 * np.float_power(x[0], -1) * np.float_power(x[1], -2) * np.float_power(intx2, -1) -1 # <= 1
  g[1] = 397.5 * np.float_power(x[0], -1) * np.float_power(x[1], -2) * np.float_power(intx2, -2) -1 # <= 1
  g[2] = 1.93 * np.float_power(x[1], -1) * np.float_power(intx2, -1) * np.power(x[3], 3) * np.float_power(x[5], -4) -1 # <= 1

  g[3] = 1.93 * np.float_power(x[1], -1) * np.float_power(intx2, -1) * np.power(x[3], 3) * np.float_power(x[5], -4) -1 # <= 1

  g[4] = (1 / (0.1 * np.power(x[5],3))) * np.power(np.power(((745*x[3])/(x[1]*intx2)),2) + 16.9e6, 0.5) - 1100

  # g[4] = (np.sqrt(np.power((745*x[3]) / (x[1]*x[2]),2)+(16.9e6)) / (110.*np.float_power(x[5], 3)) ) - 1
  
  g[5] = (1/(0.1 * np.power(x[6], 3))) * np.power((np.power(((745 * x[4]) / (x[1] * intx2)),2) + 157.5e6 ), 0.5) -850 # <=850
  g[6] = x[1] * intx2 - 40 # <= 40
  g[7] = 5 - (x[0] / x[1]) # >= 5
  g[8] = x[0] / x[1] -12 # <= 12
  g[9] = (1.5 * x[5] + 1.9) * np.float_power(x[3], -1) -1 # <= 1
  g[10] = (1.1 * x[6] + 1.9) * np.float_power(x[4], -1) -1 # <= 1
  
  return objectiveFunction, g, h

def weldedBeam(x, objFunc, g, h): # The welded beam design | 320k maxEval
  # h -> x[0] | l -> x[1] | t -> x[2] | b -> x[3]
  alpha = np.sqrt(0.25 * (np.power(x[1], 2) + np.power((x[0] + x[2]), 2))) # OK
  tau1 = (6000) / (np.sqrt(2) * x[0] * x[1]) # OK
  tau2 = ((6000 * (14 + 0.5 * x[1]) * alpha) / (2 * (0.707 * x[0] * x[1] * (np.power(x[1], 2) / 12 + 0.25 * np.power((x[0] + x[2]),2))))) # OK
  pc = 64746.022 * (1 - 0.0282346 * x[2]) * x[2] * np.power(x[3], 3) # OK

  objectiveFunction = 1.10471 * np.power(x[0], 2) * x[1] + 0.04811 * x[2] * x[3] * (14.0 + x[1]) # Cost
  g[0] = - (13600 - np.sqrt(np.power(tau1, 2) + np.power(tau2, 2) + x[1] * tau1 * tau2 / alpha)) # >= 0 # OK
  g[1] = - (30000 - 504000/(np.power(x[2], 2) * x[3])) # >= 0 OK
  g[2] = - (x[3] - x[0]) # >= 0
  g[3] = - (pc - 6000) # >= 0
  g[4] = - (0.25 - 2.1952/(np.power(x[2],3) * x[3])) # >= 0

  return objectiveFunction, g, h

def weldedBeamModObjFun1_7(x, objFunc, g, h): # The welded beam design | 320k maxEval
  '''
    P = 6000
    L = 14
    E = 30 * 10⁶
    G = 12 * 10⁶
    tau_max = 13600
    sigma_max = 30000
    delta_max = 0.25
  '''


  # h -> x[0] | l -> x[1] | t -> x[2] | b -> x[3]

  # Constants
  p = 6000
  l = 14
  e = 30 * (np.power(10, 6))
  G = 12 * (np.power(10, 6))
  tauMax = 13600
  sigmaMax = 30000
  deltaMax = 0.25

  # Defs

  m = p * (l + (x[1] / 2))
  r = np.sqrt( (np.power(x[1], 2) / 4) + np.power(((x[0] + x[2]) / 2), 2) )
  j = 2 * ( np.sqrt(2) * x[0] * x[1] * ((np.power(x[1], 2) / 12) + np.power( (x[0] + x[2] / 2), 2)) )

  tauLinha = p / (np.sqrt(2) * x[0] * x[1])
  tauDuasLinhas = (m * r) / j
  sigmaX = (6 * p * l) / ( x[3] * np.power(x[2], 2))
  deltaX = (4 * p * np.power(l, 3)) / (e * np.power(x[2], 3) * x[3])
  pcX = ((4.013 * e * np.sqrt((np.power(x[2], 2) * np.power(x[3], 6)) / 36 )) / np.power(l,2)) * (1 - (x[2] / (2 * l)) * np.sqrt(e / (4*G)))
  tauX = np.sqrt( np.power(tauLinha, 2) + 2 * tauLinha * tauDuasLinhas * (x[1] / (2*r)) + np.power(tauDuasLinhas, 2))

  
  objectiveFunction = 1.10471 * np.power(x[0], 2) * x[1] + 0.04811 * x[2] * x[3] * (14.0 + x[1])

  g[0] = tauX - tauMax # <= 0
  g[1] = sigmaX - sigmaMax # <= 0
  g[2] = x[0] - x[3] # <= 0
  g[3] = 0.10471 * np.power(x[0], 2) + 0.04811 * x[2] * x[3] * (14.0 + x[1]) - 5 # <= 0
  g[4] = 0.125 - x[0] # <= 0
  g[5] = deltaX - deltaMax # <= 0
  g[6] = p - pcX # <= 0

  return objectiveFunction, g, h

def pressureVesel(x, objFunc, g, h): # The pressure vessel design | 80k maxEval
  # Ts -> x[0] | Th -> x[1] | R -> x[2] | L -> x[3]
  # x = [0.8125, 0.4375, 42.0973, 176.6509]
  discreteSet = [0.0625, 0.125 , 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5,
                  0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.,
                  1.0625, 1.125, 1.1875, 1.25, 1.3125, 1.375, 1.4375, 1.5,
                  1.5625, 1.625, 1.6875, 1.75, 1.8125, 1.875, 1.9375, 2.,
                  2.0625, 2.125, 2.1875, 2.25, 2.3125, 2.375, 2.4375, 2.5,
                  2.5625, 2.625, 2.6875, 2.75, 2.8125, 2.875, 2.9375, 3.,
                  3.0625, 3.125, 3.1875, 3.25, 3.3125, 3.375, 3.4375, 3.5,
                  3.5625, 3.625, 3.6875, 3.75, 3.8125, 3.875, 3.9375, 4.,
                  4.0625, 4.125, 4.1875, 4.25, 4.3125, 4.375, 4.4375, 4.5,
                  4.5625, 4.625, 4.6875, 4.75, 4.8125, 4.875, 4.9375, 5.]

  x0DiscreteValue = find_nearest(discreteSet, x[0])
  x1DiscreteValue = find_nearest(discreteSet, x[1])
  objectiveFunction = 0.6224 * x0DiscreteValue * x[2] * x[3] + 1.7781 * x1DiscreteValue * np.power(x[2], 2) + 3.1661 * np.power(x0DiscreteValue, 2) * x[3] + 19.84 * np.power(x0DiscreteValue, 2) * x[2] # Weight
  g[0] = - (x0DiscreteValue - 0.0193 * x[2]) # >= 0
  g[1] = - (x1DiscreteValue - 0.00954 * x[2]) # >= 0
  g[2] = - (np.pi * np.power(x[2], 2) * x[3] + 4/3 * np.pi * np.power(x[2], 3) - 1296000) # >= 0
  g[3] = - (-x[3] + 240) # >= 0

  return objectiveFunction, g, h

def cantileverBeam(x, objFunc, g, h): # The cantilever beam design | 35k maxEval
  # Bi -> x[0] | Hi -> x[1]
  # x = [3, 60, 3.1, 55, 2.6, 50, 2.2837, 45.5507, 1.7532, 35.0631] # Best solution Bernardino 08
  bDiscreteSet = [2.4, 2.6, 2.8, 3.1]
  hDiscreteSet = [45., 50., 55., 60.]
  # x[0] => b1 # Integer
  # x[1] => h1 # Integer

  # x[2] => b2 # Discrete
  # x[3] => h2 # Discrete
  # x[4] => b3 # Discrete
  # x[5] => h3 # Discrete
  # x[6] => b4 # Continuous
  # x[7] => h4 # Continuous
  # x[8] => b5 # Continuous
  # x[9] => h5 # Continuous

  intx0 = np.round(x[0])
  intx1 = np.round(x[1])

  discretex2 = find_nearest(bDiscreteSet, x[2])
  discretex3 = find_nearest(hDiscreteSet, x[3])
  discretex4 = find_nearest(bDiscreteSet, x[4])
  discretex5 = find_nearest(hDiscreteSet, x[5])


  # objectiveFunction = 100*(x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9]) # Volume
  objectiveFunction = 100*(intx0*intx1 + discretex2*discretex3 + discretex4*discretex5 + x[6]*x[7] + x[8]*x[9]) # Volume

  g[0] = 6*25000000/(intx0*(intx1*intx1)) - 14000
  g[1] = 6*20000000/(discretex2*(discretex3*discretex3)) - 14000
  g[2] = 6*15000000/(discretex4*(discretex5*discretex5)) - 14000
  g[3] = 6*10000000/(x[6]*(x[7]*x[7])) - 14000
  g[4] = 6*5000000/(x[8]*(x[9]*x[9])) - 14000
  g[5] = (intx1/intx0) - 20
  g[6] = (discretex3/discretex2) - 20
  g[7] = (discretex5/discretex4) - 20
  g[8] = (x[7]/x[6]) - 20
  g[9] = (x[9]/x[8]) - 20

  aj1 = intx0*(np.power(intx1, 3))/12
  aj2 = discretex2*(np.power(discretex3, 3))/12
  aj3 = discretex4*(np.power(discretex5, 3))/12
  aj4 = x[6]*(np.power(x[7], 3))/12
  aj5 = x[8]*(np.power(x[9], 3))/12
  g[10] = (0.0025/3)*( (1000000/aj5) + (7000000/aj4) + (19000000/aj3) + (37000000/aj2) + (61000000/aj1) ) - 2.7

  # print("Objective Function: {}".format(objectiveFunction))
  # print("g: {}".format(g))
  # sys.exit("Vlw")
  return objectiveFunction, g, h

def executeFunction(function, x, objFunc, g, h):
  if function == 21:
    objFunc, g, h = tensionCompressionSpring(x, objFunc, g, h) # 36k maxEval
  elif function == 22 or function == 222:
    objFunc, g, h = speedReducer(x, objFunc, g, h) # 36k maxEval
  elif function == 23:
    objFunc, g, h = weldedBeam(x, objFunc, g, h) # 320k maxEval
  elif function == 233:
    objFunc, g, h = weldedBeamModObjFun1_7(x, objFunc, g, h) # 320k maxEval
  elif function == 24:
    objFunc, g, h = pressureVesel(x, objFunc, g, h) # 80k maxEval
  elif function == 25:
    objFunc, g, h = cantileverBeam(x, objFunc, g, h) # 35k maxEval
  else:
    sys.exit("Function not defined.")
  return objFunc, g, h

def tam(x):
  return 4 + np.floor(3 * np.log(x))
# if __name__ == '__main__':
#   executeFunction()




#  Comparações com paper DeMelo

# Welded beam 20k
# Speed Reducer 30k
# Three bar truss 10k
# Pressure vessel 30k
# Tension/compression spring 20k

# Function 21 => Tension compression spring | 20k
# Function 222 => Speed reducer | 30k
# Function 233 => Welded beam | 20k
# Function 24 => Pressure vessel | 30k
# Function 25 => Cantilever beam


# Shen, Zhu, Niu & Wu 2009
# Ray & Liew 2003
# Coelho 2010



# *---*---*---*---*---*---*---*---*---* Pandas Personal Helper ---*---*---*---*---*---*---*---*---*---*---*

# g(x) <= V -> V - g(x) <= 0
# g(x) >= V -> g(x) - V >= 0 

# def rosen():
# 	if function == 91:  # RosenbrockFunction
# 	sumRosen = 0
# 	for k in range(nSize - 1):
# 		sumRosen = sumRosen + 100*np.power((self.individuals[i].n[k+1] - np.power(self.individuals[i].n[k], 2)), 2) + np.power((1 - self.individuals[i].n[k]), 2)
# 	self.individuals[i].objectiveFunction[0] = sumRosen
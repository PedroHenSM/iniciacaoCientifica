import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
from scipy.stats import ranksums # Wilcoxon signed rank test
import os
import sys
from pathlib import Path

# PROBLEMS_TYPE = "Classical Engineering" # Trusses | Classical Engineering | All
PROBLEMS_TYPE = "Trusses" # Trusses | Classical Engineering | All
PP_MEASUREMENT = "Best" # Best | Median | Average | St. Dev | Worst
RUN_WILCOXON_TEST = False

class Analysis(object):
  def __init__(self, columns=None, rows=None, currentFunction=0, functionsName=None, functionsMaxEval=None, problemsSize=None):
    # self.table = [[ 0 for _ in range(len(columns))] for _ in range(len(rows))] if columns is not None else None
    self.table = []
    self.columns = columns if columns is not None else []
    self.rows= rows if rows is not None else []
    self.functionsName = functionsName if functionsName is not None else []
    self.functionsMaxEval = functionsMaxEval if functionsMaxEval is not None else {}
    self.currentFunction = currentFunction
    self.problemsSize= problemsSize if problemsSize is not None else 0

  # Makes print(individual) a string, not a reference (memory adress)
  def __repr__(self):
    return str(self.__dict__) + "\n"

def getNameIndividualToChoose(individual):
  individualToChoose = None
  if individual == "Hof":
    individualToChoose = "Hall of fame"
  elif individual == "Last":
    individualToChoose = "Last individual"
  elif individual == "Last factible":
    individualToChoose = "Last factible individual"
  else:
    sys.exit("Individual to pick not defined.")
  return individualToChoose

def readResults():
  # algorithms = ["DE", "CMAES"]
  algorithms = ["CMAES"]
  if PROBLEMS_TYPE == "Trusses":
    functions = [110, 125, 160, 172]
    # functions = [110]
    # functions = [110, 160]
    # functions = [160]
  elif PROBLEMS_TYPE == "Classical Engineering":
    functions = [21, 22, 23, 24, 25] # Classical engineering problems
    # functions = [21] # Classical engineering problems
  elif PROBLEMS_TYPE == "All":
    # functions = [21, 22, 23, 24, 25, 110, 125, 160, 172, 1942] # All problems
    functions = [21, 22, 23, 24, 25, 110, 125, 160, 172] # All problems

  # weights = ["Linear", "Superlinear", "Equal"] # CMA-ES weights parameters
  weights = ["Superlinear"] # CMA-ES weights parameters

  # individualToChoose = ["Last factible", "Hof"] # Types of individual to choose for analysis
  individualToChoose = ["Hof"] # Types of individual to choose for analysis
  trussCases = ["c", "d"] # Continuous or Discrete
  functionsName = [] # Saves functions name a strings ..:  ['110c', '110d', '125d']
  # constraintHandlingMethods = ["APM"]
  constraintHandlingMethods = ["DEB", "APM"]



  seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  , 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
  # seeds = [1, 2]

  basePath = Path(__file__).parent
  solutions = []

  restartCriterias = {
    'tolFun': 0,
    'tolX': 0,
    'noEffectAxis': 0,
    'noEffectCoord': 0,
    'conditionCov': 0,
  }

  # Dictionary containing functions to be analysed with max evaluations allowed
  functionsMaxEval = {
    # Trusses
    110: {
      'c': 280000,
      'd': 90000
      # 'c': False,
      # 'd': 240000
    },
    125: {
      # 'c': 240000,
      'c': False,
      'd': 20000
    },
    160: {
      'c': 12000,
      # 'c': 800000,
      'd': False
    },
    172: {
      'c': 35000,
      'd': False
    },
    1942: {
      'c': False,
      'd': False
    },
    # Classical engineering problems
    21: 36000,
    22: 36000,
    23: 320000,
    24: 80000,
    25: 35000
  }

  for f in functions: # Functions
    for c in trussCases:
      dataOfFunction = []
      for a in algorithms: # Algorithms
        for p in constraintHandlingMethods: # Constraint handling methods
          for keyW, w in enumerate(weights):
            for i in individualToChoose:
              for s in (seeds): # Seeds
                tempData = []
                filePath = None
                if a == "CMAES":
                  if str(f)[0] == "1": # Truss problems
                    # Get the max evaluations
                    m = functionsMaxEval[f][c]
                    if m:
                      filePath = (basePath / "../results/functions/f{}/{}_f{}_c{}_p{}_w{}_m{}_s{}.dat".format(f, a, f, c, p, w, m, s)).resolve()
                  elif str(f)[0] == "2": # Engineering problems
                    # Get the max evaluations
                    m = functionsMaxEval[f]
                    if m:
                      filePath = (basePath / "../results/functions/f{}/{}_f{}_p{}_w{}_m{}_s{}.dat".format(f, a, f, p, w, m, s)).resolve() # Original  code
                elif a == "DE" and keyW == 0: # Only make analysis once, since weights parameters doesn't change DE results
                  if str(f)[0] == "1": # Truss problems
                    # Get the max evaluations
                    m = functionsMaxEval[f][c]
                    if m:
                      filePath = (basePath / "../results/functions/f{}/{}_f{}_c{}_p{}_w{}_m{}_s{}.dat".format(f, a, f, c, p, w, m, s)).resolve() # Original  code
                  elif str(f)[0] == "2": # Engineering problems
                    # Get the max evaluations
                    m = functionsMaxEval[f]
                    if m:
                      filePath = (basePath / "../results/functions/f{}/{}_f{}_p{}_m{}_s{}.dat".format(f, a, f, p, m, s)).resolve() # Original  code
                # Found path
                if filePath is not None:
                  # Read file
                  if str(f)[0] == "1" and not str(f)+c in functionsName: # Trusses
                    functionsName.append(str(f)+c)
                  elif str(f)[0] == "2" and not str(f) in functionsName: # Engineering problems
                    functionsName.append(str(f))
                  file = open(filePath)
                  countFactibleInd = 0
                  while True:

                    buffer = file.readline()
                    # print(buffer)
                    if 'Restart causes:' in buffer:
                      buffer = buffer.split("'")
                      # print("Spplited buffer: {}".format(buffer))
                      restartCriterias[buffer[1]] +=1 
                      # sys.exit("ok")
                    if getNameIndividualToChoose(i) in buffer: # Find individual for data analysis
                      hasFactibleSolution = True
                      updateBestIndividual = False
                      buffer = file.readline() # Read one more
                      if len(buffer) > 1: # Buffer is not empty
                        buffer = buffer.split(" ")
                        # Verify if solution is feasible
                        if p == "DEB":
                          if float(buffer[1]) != 0:
                            hasFactibleSolution = False
                            # print("Individual infactible. Problem: {}: {}+{}+s{}".format(f, a, p, s))
                        elif p == "APM":
                          if float(buffer[0]) != float(buffer[1]):
                            hasFactibleSolution = False
                            # print("Individual infactible. Problem: {}: {}+{}+s{}".format(f, a, p, s))
                        
                        # Only saves individual if its factible
                        if hasFactibleSolution:
                          # print("Another factible individual: {}: {}+{}+s{}".format(f, a, p, s))
                          for key, item in enumerate(buffer):
                            # First individual to be inserted on tempData, just insert it
                            if countFactibleInd == 0:
                              tempData.append(float(item))
                            else:
                              if key == 0: # First buffer item
                                if float(item) < tempData[-len(buffer)]: # Ojbective function from new individual is better than old one
                                  # print("Found a better one. Problem: {}: {}+{}+s{}".format(f, a, p, s))
                                  # print("tempData before cleaning: {}".format(tempData))
                                  tempData = tempData[:-len(buffer)] # Remove last individual
                                  # print("tempData after cleaning: {}".format(tempData))
                                  updateBestIndividual = True
                                else: # If new individual is worst than old one, dont need go through the buffer array
                                  break
                              if updateBestIndividual: # If individual has to be updated
                                # print("Found a better one. Problem: {}: {}+{}+s{}".format(f, a, p, s))
                                tempData.append(float(item))
                          countFactibleInd += 1

                    # Read time and go to next one
                    elif "CPU time used" in buffer:
                      # Only get CPU time and algorithm name that algorithm got at least one factible solution
                      if countFactibleInd != 0: 
                        buffer = file.readline()
                        tempData.append(float(buffer))
                        if a == "CMAES":
                          if str(f)[0] == "1": # Truss problems
                            tempData.append("Problem{}{}_{} {} {} + {} - m{}".format(f, c, a, w, i, p, m))
                          elif str(f)[0] == "2": # Engineering problems
                            tempData.append("Problem{}_{} {} {} + {} - m{}".format(f, a, w, i, p, m))
                        elif a == "DE":
                          if str(f)[0] == "1": # Truss problems
                            tempData.append("Problem{}{}_{} {} + {} - m{}".format(f, c, a, i, p, m))
                          elif str(f)[0] == "2": # Engineering problems
                            tempData.append("Problem{}_{} {} + {} - m{}".format(f, a, i, p, m))
                        
                        dataOfFunction.append(tempData)
                      break
      # Each list from solutions contains the data for each function
      if isTruss(f):
        # Only appends if has results (dataOfFunctions is not empty)
        if dataOfFunction:
          solutions.append(dataOfFunction)
    if not isTruss(f):
      # Only appends if has results (dataOfFunctions is not empty)
      if dataOfFunction:
        solutions.append(dataOfFunction)

  # print("Functions name: {}".format(functionsName))
  # sys.exit("ok")

  #   restartCriteriasPercentage = {
  #   'tolFun': 0,
  #   'tolX': 0,
  #   'noEffectAxis': 0,
  #   'noEffectCoord': 0,
  #   'conditionCov': 0,
  # }

  # totalRestarts = sum(restartCriterias.values())
  # for criteria, value in restartCriterias.items():
  #   restartCriteriasPercentage[criteria] = np.round((value / totalRestarts) * 100, 2)
  # print("Valores absolutos: {} \\\\".format(restartCriterias))
  # print("Porcentagem: {} \\\\".format(restartCriteriasPercentage))
  # sys.exit("Exiting")
  return solutions, functionsName, functionsMaxEval

def getExtraResultsNew(np2dArr, rowsTitles, function):
  # Trusses problems
  if function == "110c": # 10 bar truss continuous
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 5060.87, 5060.92, 5061.98, 3.93e+00, 5076.70, float("NaN")] 
    duvde = ["DUVDE - m280000", 5060.85, float("NaN"), 5067.18, 7.94e+00, 5076.66, float("NaN")] 
    apm = ["APM*", 5069.08, float("NaN"), 5091.43, float("NaN"), 5117.39, float("NaN")]
    # Tese
    # apmErica = ["PSO + APM Erica - m280000", 5060.9770, 5076.8523, 5072.2409, 4.8201e+01, 5087.7249, float("NaN")]
    # Executados por mim
    apmErica = ["PSO + APM - m280000", 5060.947307, 5069.038276, 5069.816972, 4.629297e+01, 5087.805227, float("NaN")]
    apmRafael = ["GA + APM - m280000", 5061.935, float("NaN"), 5069.193, 6.93e+00, 5086.705, float("NaN")]

    # Append name of algorihtms
    # rowsTitles.append(smde17[0])
    rowsTitles.append(duvde[0])
    # rowsTitles.append(apm[0])
    rowsTitles.append(apmErica[0])
    rowsTitles.append(apmRafael[0])

    # Append reults on np2darr
    # np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
    np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
    # np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmErica[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmRafael[1:]], axis=0)
  elif function == "110d": # 10 bar truss discrete
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 5490.74, 5490.74, 5495.99, 1.13e+01, 5529.30, float("NaN")]
    duvde = ["DUVDE - m90000/24000", 5562.35, float("NaN"), 5564.90, 0.6, 5565.04, float("NaN")]
    apm = ["APM*", 5490.74, float("NaN"), 5545.48, float("NaN"), 5567.84, float("NaN")]
    # Tese
    # apmErica = ["PSO + APM Erica - m90000", 5509.7173, 5528.0869, 5628.8689, 1.4468e+03, 6593.1205, float("NaN")]
    # Executados por mim
    # apmErica = ["PSO + APM - m180000", 5509.7173, 5521.244867, 5594.236423, 8.671090e+02, 6347.034801, float("NaN")]
    apmErica = ["PSO + APM - m90000", 5509.717373, 5583.433879, 5608.054915, 9.398600e+02, 6463.898099, float("NaN")]
    apmRafael = ["GA + APM - m90000", 5538.321, float("NaN"), 5543.38, 8.41e+00, 5559.598, float("NaN")]

    # Append name of algorihtms
    # rowsTitles.append(smde17[0])
    rowsTitles.append(duvde[0])
    # rowsTitles.append(apm[0])
    rowsTitles.append(apmErica[0])
    rowsTitles.append(apmRafael[0])

    # Append reults on np2darr
    # np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
    np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
    # np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmErica[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmRafael[1:]], axis=0)
  elif function == "125c": # 25 bar truss continuous
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 484.06, 484.07, 484.07, 0.0107, 484.10, float("NaN")]
    duvde = ["DUVDE - m240000", 484.05, float("NaN"), 484.05, 1e-05, 484.05, float("NaN")]
    
    # Append name of algorihtms
    # rowsTitles.append(smde17[0])
    rowsTitles.append(duvde[0])

    # Append reults on np2darr
    # np2dArr = np.append(np2dArr, [smde17[1:]], axis=0) # Appends values of given list (remove first index, which is the string)
    np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
  elif function == "125d": # 25 bar truss discrete
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 484.85, 485.05, 485.44, 0.693, 487.13, float("NaN")]
    duvde = ["DUVDE - m20000", 485.90, float("NaN"), 498.44, 7.66e+00, 507.77, float("NaN")]
    apm = ["APM*", 485.85, float("NaN"), 485.97, float("NaN"), 490.74, float("NaN")]
    # Tese
    # apmErica= ["PSO + APM Erica - m?", 484.8541, 485.0487, 485.97338, 1.3601e+01, 496.2520, float("NaN")]
    # Executados por mim
    apmErica= ["PSO + APM - m20000", 484.854179, 485.048797, 486.289066, 2.526693e+01, 510.336915, float("NaN")]
    apmRafael= ["GA + APM - m20000", 484.8542, float("NaN"), 485.214, 5.89e-01, 487.3456, float("NaN")]

    # Append name of algorihtms
    # rowsTitles.append(smde17[0])
    rowsTitles.append(duvde[0])
    # rowsTitles.append(apm[0])
    rowsTitles.append(apmErica[0])
    rowsTitles.append(apmRafael[0])

    # Append reults on np2darr
    # np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
    np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
    # np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmErica[1:]], axis=0) # Appends values of given list (remove first index, which is the string)
    np2dArr = np.append(np2dArr, [apmRafael[1:]], axis=0)
  elif function == "160c": # 60 bar truss continuous
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 308.94, 309.42, 309.49, 0.464, 311.21, float("NaN")]
    duvde = ["DUVDE - m150000", 309.44, float("NaN"), 311.54, 1.46e+00, 314.70, float("NaN")]
    apm = ["APM", 311.87, float("NaN"), 333.01, float("NaN"), 384.19, float("NaN")]
    # Tese
    # apmErica = ["PSO + APM Erica - m?", 291.2477, 302.0956, 319.2243, 1.9759e+02, 436.8457, float("NaN")]
    # Executados por mim
    apmErica = ["PSO + APM - m12000", 311.717441, 319.938600, 327.910086, 9.568748e+01, 380.653120, float("NaN")]
    apmRafael = ["GA + APM - m12000", 308.6482 , float("NaN"), 330.633 , 1.38e+01, 360.4936, float("NaN")]

    # Append name of algorihtms
    # rowsTitles.append(smde17[0])
    rowsTitles.append(duvde[0])
    # rowsTitles.append(apm[0])
    rowsTitles.append(apmErica[0])
    rowsTitles.append(apmRafael[0])

    # Append reults on np2darr
    # np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
    np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
    # np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmErica[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmRafael[1:]], axis=0)
  elif function == "160d": # 60 bar truss discrete
    # Define algorithms and literature results
    smde = ["SMDE k=2* (2017)", 312.73, 314.20, 315.12, 3.98e+00, 335.88, float("NaN")]

    # Append name of algorihtms
    rowsTitles.append(smde[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
  elif function == "172c": # 72 bar truss continuous
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 379.62, 379.63, 379.65, 0.0341, 379.73, float("NaN")]
    duvde = ["DUVDE - m35000", 379.66, float("NaN"), 380.42, 0.572, 381.37, float("NaN")]
    apm = ["APM", 387.04, float("NaN"), 402.59, float("NaN"), 432.95, float("NaN")]
    # Tese
    # apmErica = ["PSO + APM Erica - m?", 379.65373, 379.74019, 383.87949, 1.0493e+02, 475.87770, float("NaN")]
    # Executador por mim
    apmErica = ["PSO + APM - m35000", 379.652909, 379.731173, 379.745805, 5.449887e-01, 380.204020, float("NaN")]
    apmRafael = ["GA + APM - m35000", 383.0324, float("NaN"), 388.537, 3.26e+00, 394.8527, float("NaN")]

    # Append name of algorihtms
    # rowsTitles.append(smde17[0])
    rowsTitles.append(duvde[0])
    # rowsTitles.append(apm[0])
    rowsTitles.append(apmErica[0])
    rowsTitles.append(apmRafael[0])

    # Append reults on np2darr
    # np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
    np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
    # np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmErica[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmRafael[1:]], axis=0)
  elif function == "172d": # 72 bar truss discrete
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 385.54, 386.81, 386.91, 1.05e+00, 389.21, float("NaN")]

    # Append name of algorihtms
    rowsTitles.append(smde17[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
  elif function == "1942c": # 942 bar truss continuous
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 149932.00, 171218.50, 174369.63, 1.72e+04, 230139.00, float("NaN")]

    # Append name of algorihtms
    rowsTitles.append(smde17[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
  elif function == "1942d": # 942 bar truss discrete
    # Define algorithms and literature results
    smde17 = ["SMDE k=2* (2017)", 153010.00, 181899.00, 181127.23, 1.73e+04, 216514.00, float("NaN")]

    # Append name of algorihtms
    rowsTitles.append(smde17[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [smde17[1:]], axis=0)
  # Other Engineering Problems
  elif function == "21": # Tension/compression spring
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 0.012666, 0.012892, 0.013131, 6.28e-4, 0.015318, 50] 
    apmSpor = ["GA + APM SPOR - m36000", 0.012667602164, float('NaN'), 0.013748492439, float('NaN'), 0.017093902154, float('NaN')]
    # Tese
    # apmMed3 = ["PSO + APM MED 3 - m36000", 0.01266, 0.01312, 0.01389 , 9.1731e-03, 0.01777, 35]
    # Executados por mim
    apmMed3 = ["PSO + APM MED 3 - m36000", 0.012666, 0.012895, 0.013597 , 7.200257e-03, 0.017406, 30]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == "22": # Speed reducer
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 2996.3483, 2996.3495, 2996.3501, 7.45e-3, 2996.3599, 50] 
    apmSpor = ["GA + APM SPOR - m36000", 2996.34850933205, float('NaN'), 2996.35243640334, float('NaN'), 2996.36609677358, float('NaN')]
    # Tese
    # apmMed3 = ["PSO + APM MED 3 - m36000", 2996.3622 , 2996.3780, 2999.6083, 3.4911e+01, 3016.7808, 35]
    # Executados por mim
    apmMed3 = ["PSO + APM MED 3 - m36000", 2996.356314 , 2996.381070, 3000.659737, 3.429298e+01, 3016.777031, 29]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == "23": # Welded beam
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 2.38335, 2.92121, 2.99298, 2.02e-1, 4.05600, 50] 
    apmSpor = ["GA + APM SPOR - m?", 2.38113481849464 , float('NaN'), 2.58228221674671 , float('NaN'), 3.20898593483156, float('NaN')]
    # Tese
    # apmMed3 = ["PSO + APM MED 3 - m320000", 2.38114 , 2.43315, 2.67102, 2.0656e+00, 3.46638, 35] #nRun: 32k
    # Executados por mim
    # apmMed3 = ["PSO + APM MED 3 - m32000", 2.381162 , 2.521305, 2.763862, 2.663334e+00, 4.127299, 30] #nRun: 32k
    # Executados por mim
    apmMed3 = ["PSO + APM MED 3 - m320000", 2.381146 , 2.512969, 2.663524, 1.984252e+00, 3.585389, 30] #nRun:320k

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == "24": # Pressure Vesel
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 6059.855, 6426.710, 6545.126, 1.24E+2, 7388.160, 50] 
    apmSpor = ["GA + APM SPOR - m80000", 6059.73045731256 , float('NaN'), 6581.18398763114 , float('NaN'), 7333.93495942434, float('NaN')]
    # Tese
    # apmMed3 = ["PSO + APM MED 3 - m80000", 6059.7143 , 6370.7797, 6427.6676, 2.6221e+03, 7544.4925, 35]
    # Executados por mim
    apmMed3 = ["PSO + APM MED 3 - m80000", 6059.714360, 6370.779797, 6504.021887, 2.701649e+03, 7544.492518, 30]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == "25": # Cantilever beam
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 64834.70, 74987.16, 76004.24, 6.93E+3, 102981.06, 50] 
    # best median avg std worst
    apmSpor = ["GA + APM SPOR - m35000", 64599.980343 , float('NaN'), 66684.276327 , float('NaN'), 72876.210779, float('NaN')]
    # Tese
    # apmMed3 = ["PSO + APM MED 3 - m35000", 64578.271, 68294.702, 71817.816, 1.0431e+05, 173520.325, 35]
    # Executados por mim
    apmMed3 = ["PSO + APM MED 3 - m35000", 64578.225361, 68943.497363, 68873.965835, 1.703295e+04, 76243.501901, 30]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  else:
    sys.exit("Function not defineda.")
  return np2dArr, rowsTitles

def makeAnalysisNew(solutions, functions, analysis):
  # Read solutions and store on dataframe
  # case="discrete"
  # print("solutions len: {}".format(len(solutions)))
  # print("solutions[idx]: {}".format(solutions[0]))
  # print("functions:  {}".format(functions))
  # print("solutions: {}".format(solutions))
  # sys.exit('thx')
  for key, problem in enumerate(solutions):
    # Verifies if solution array is not empty
    if problem:
      # print("key: {}".format(key))
      df = pd.DataFrame.from_records(problem)
      titles = ["Objective Function", "ViolationSum/Fitness"]
      # Define name for columns of the project variables
      for i in range(len(df.columns)-4):
        titles.append("DesignVariables{}".format(i))

      # Define name for the last 2 columns
      titles.append("CPU Time(s)")
      titles.append("Algorithm")

      # Name columns
      df.columns = titles

      # Drop columns that contains word | Could be done with df.filter(regex)
      # df.drop([col for col in df.columns if "DesignVariables" in col], axis=1, inplace=True)
      # df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
      df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)



      # # OBJFUNC T60C: 309.115 violationSum: 309.115
      # bestT60C = [1.14806, 2.0408, 0.501062, 1.85304, 1.75507, 0.585083, 1.9266, 1.8511, 1.00357, 1.88154, 1.86348, 0.5, 2.02915, 1.24715, 1.09826, 0.637386, 0.730659, 0.972064, 1.12724, 1.14845, 1.05482, 1.05667, 0.654224, 1.02509, 1.25693]

      # print("df:\n{}".format(df))
      # bestObj = df[df['Objective Function'] == df['Objective Function'].min()]
      # bestObj1 = df.iloc[df['Objective Function'].idxmin()]
      # print(bestObj)
      # print("\n\n")
      # print(bestObj1)
      # print("\n\n")
      
      # sys.exit()
      # print("df min ObjFunc:\n{}".format(df[df['Objective Function'] == df['Objective Function'].min()]))
      # sys.exit('okz')
      # # Group by algorithm and calcualtes statistical measures
      grouped = df.groupby(["Algorithm"])
      # print("grouped:\n{}".format(grouped))
      # sys.exit('ok')

      # Imprime as variáveis de projeto de cada problema, por algoritmo
      pd.set_option('display.max_columns', None)
      pd.set_option('display.max_rows', None)
      pd.set_option('display.width', None)
      pd.set_option('display.max_colwidth', None)
      # groupedMin = grouped.min()
      # print(groupedMin)

      # print('Imprime Grouped MIN (FIXED)')
      print(problem[key][-1])
      groupedMin = df.loc[df.groupby(['Algorithm'])['Objective Function'].idxmin()]
      print(groupedMin)
      print()


      # sys.exit('ok')
      # Bloco abaixo funciona
      myDict = {}
      algorithmArrayDf = grouped['Objective Function'].apply(list)
      for idx, item in enumerate(algorithmArrayDf.index):
        # print('aaaaaa {}'.format(idx))
        # print('aaaaaa {}'.format(item))
        myDict[item] = algorithmArrayDf.values[idx]

      # # Print abaixo funcional!
      # print("my dict:\n{}".format(myDict))
      # print("\nObjective functions for each seed\n")
      # for keys, values in myDict.items():
      #   print(keys)
      #   print(values)
      # # sys.exit('ok')
      # # Bloco acima funciona
      

      # grouped = df.groupby(["Algorithm"])
      grouped = grouped.agg(["min", "median", "mean", "std", "max", "size"]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works
      # print(grouped)
      # Generates 2d np array from pandas grouped df
      np2dArr = grouped.to_numpy()
      # Get row and columns titles
      rowsTitles = list(grouped.index.values) 
      columnsTitles = list(grouped.columns.levels[1])
      # Get extra results (appending other literature algorithms on 2d np arr)
      # np2dArr, rowsTitles = getExtraResultsNew(np2dArr, rowsTitles, functions[key])
      # Generates .text table
      # if functions[key] == '110c' or functions[key] == '110d' or functions[key] == '125d' or functions[key] == '160c' or functions[key] == '172c':
      # print("analysis: {}".format(analysis))
      # sys.exit('obrigado')
      # Imprime tabelas com os resultados
      # np2dArrToLatexNew(np2dArr, columnsTitles, rowsTitles, analysis)
  # sys.exit('ok')

def isTruss(function):
  if str(function)[0] == "1": # Truss problem
    return True
  return False

def getItemIndex(item, analysisTable):
  indexOf = None
  try:
    indexOf = analysisTable.rows.index(item)
  except ValueError:
    indexOf = None;
  return indexOf

# Removes problem name from rows titles
def adjustRowTitle(rowsTitles): 
  # Problem110d_CMAES Linear Hof + APM
  # Removes underscore from algorithm names and maxEvaluations after '-'
  algorithmsMaxEval = {}
  for i in range(len(rowsTitles)):
    splitted = rowsTitles[i].split(" - ")
    if len(splitted) == 1:
      title = splitted[0] # Algorithm name
    elif len(splitted) == 2:
      title = splitted[0] # Algorithm name
      maxEval = splitted[1][1:] # Max evaluations for this problem
      if 'Problem' in title: # My algorithms executed by myself
        algorithmsMaxEval['Max evaluations'] = maxEval
      else: # OAlgorithms from literature
        algorithmsMaxEval[title] = maxEval
    splitted = title.split("_")
    # splitted = rowsTitles[i].split("_")
    if len(splitted) == 1:
      rowsTitles[i] = splitted[0]
    elif len(splitted) == 2:
      rowsTitles[i] = splitted[1]
    else:
      sys.exit("Unexpected splitted size. Please verify method")
  # return rowsTitles, maxEval
  return rowsTitles, algorithmsMaxEval

def np2dArrToLatexNew(np2dArr, columnsTitles, rowsTitles, analysisTable):
  # Current function being executed
  currentFunctionName = analysisTable.functionsName[analysisTable.currentFunction]

  # Current function is a truss, print if problem is continuous or discrete
  trussCase = None
  if isTruss(currentFunctionName):
    if currentFunctionName[-1] == "d": # Discrete problem
      trussCase = "Discrete"
    else: # Continuous problem
      trussCase = "Continuous"

  # Gets caption biased on string before underscore and truss case (if truss problem)
  caption = ("{} {}".format(renameProblem(rowsTitles[0].split("_")[0]), trussCase)) if trussCase is not None else ("{}".format(renameProblem(rowsTitles[0].split("_")[0]))) 

  # Rename column titles
  columnsTitles = renameColumnsTitles(columnsTitles)
  rowsTitles, algorithmsMaxEval = adjustRowTitle(rowsTitles)

  # Generate files for performance profiles
  ppProblemsFilePath = "../results/analysis/performanceProfiles/ppProblems.dat"
  ppProblemsFile = open(ppProblemsFilePath, "a")
  # Running the first problem, clear ppProblems.dat file 
  if not analysisTable.currentFunction: 
    ppProblemsFile.truncate(0)
  # Write current problem name
  ppProblemsFile.write('{}\n'.format(caption.replace(" ", "_")))
  ppProblemsFile.close()

  ppAlgorithmsFilePath = "../results/analysis/performanceProfiles/ppAlgorithms.dat"
  ppAlgorithmsFile = open(ppAlgorithmsFilePath, "w")
  for title in rowsTitles:
    ppAlgorithmsFile.write('{}\n'.format(title))
  ppAlgorithmsFile.close()

  for title in rowsTitles:
    itemIdx = getItemIndex(title, analysisTable)
    if itemIdx is None: # Item doesn't exist on analysis.table.rows, must insert it
      dataRow = [ 0 for _ in range(len(analysisTable.columns))] # Creates an array with zeros
      analysisTable.rows.append(title)
      analysisTable.table.append(dataRow)
      # print("title: {}".format(title))
      # print("dataRow: {}".format(dataRow))
    # Now exists, increments run count
    itemIdx = getItemIndex(title, analysisTable)
    analysisTable.table[itemIdx][-1] += 1
    # else: # Item already exists. Increment run count
    #   analysisTable.table[itemIdx][-1] += 1
  # sys.exit("ok")
  
  # List containing minimum value of each column
  minOfEachColumn = np.nanmin(np2dArr, axis=0)

  # Define print format, table aligmnent and round rpecision
  precRound = 4 # Precision round parameter for highlightning minimum value from column
  if isTruss(currentFunctionName): # Truss problems, use only 2 floating points
    precRound = 2
  precFormat = "e" # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  if (type(precRound) == int):
    precFormat = ".{}f".format(precRound) # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  
  texTableAlign = "r " # Align on right (could be l, c, r and so on)
  tabAlign = "{{@{{}} l | {} @{{}}}}".format(texTableAlign*len(columnsTitles)) 

  # Begin of table structure .tex
  print("\\begin{table}[h]")
  print("\\centering")
  print("\\resizebox{0.8\\textwidth}{!}{")
  print("\\begin{minipage}{\\textwidth}")
  print("\\caption{{{}}}".format(caption))
  
  # print("\\vspace{{{}}}".format("0.5cm"))
  print("\\begin{{tabular}} {}".format(tabAlign))
  print("\\hline")
  print("&", end=" ")
  print(*columnsTitles, sep=" & ", end=" \\\ ") # Printing columnTitles
  print("\n\\hline")


  for rowIndex, row in enumerate(np2dArr):
    # Opens/creates file for each algorithm
    ppAlgorithmFilePath = ("../results/analysis/performanceProfiles/{}.dat".format(rowsTitles[rowIndex]))
    ppAlgorithmFile = open(ppAlgorithmFilePath, "a")
    # Running the first problem, clear ppProblems.dat file 
    if not analysisTable.currentFunction: 
      ppAlgorithmFile.truncate(0)
    
    # Printing row titles
    print(rowsTitles[rowIndex], end=" & ")

    # print("********************* rowTitles: {} ********************".format(rowsTitles))
    for columnIndex, item in enumerate(row):
      # print("columnsTitles[columnIndex]: {}".format(columnsTitles[columnIndex]))
      roundedValue = np.round(np2dArr[rowIndex][columnIndex], precRound)
      if columnsTitles[columnIndex] == PP_MEASUREMENT:
        ppAlgorithmFile.write("{}\t{}\n".format(caption.replace(" ", "_"), roundedValue ))
      # print("********************* item: {} ********************".format(item))
      # Verifies if current item from matrix is the minimum, for highlightning
      if roundedValue == np.round(minOfEachColumn[columnIndex], precRound):
        # Increments algorithm who performed better
        idxItem = getItemIndex(rowsTitles[rowIndex],analysisTable)
        analysisTable.table[idxItem][columnIndex] +=1
        # Last item of row (fr) doesn't need bold, its printed without floating points and uses line break 
        if columnIndex == len(row) -1: # TODO Use name of column instead of indexes for more readability
          print("{:{prec}}".format(item, prec=".0f"), end=" \\\ ")
        # Third to last item of row (std) its printed used scientific notation 
        elif columnIndex == len(row) -3: # TODO Use name of column instead of indexes for more readability
          print("\\textbf{{{:{prec}}}}".format(item, prec="0.2e"), end=" & ")
        else:
          print("\\textbf{{{:{prec}}}}".format(item, prec=precFormat), end=" & ")
      else:
        # Last item of row (fr) its printed without floating points and uses line break
        if columnIndex == len(row) -1: # TODO Use name of column instead of indexes for more readability
          print("{:{prec}}".format(item, prec=".0f"), end=" \\\ ")
        # Third to last item of row (std) its printed used scientific notation 
        elif columnIndex == len(row) -3: # TODO Use name of column instead of indexes for more readability
          print("{:{prec}}".format(item, prec="0.2e"), end=" & ")
        else:
          print("{:{prec}}".format(item, prec=precFormat), end=" & ")
    print()
    ppAlgorithmFile.close()

  print("\\end{tabular}")
  print("\\end{minipage}}")
  print("\\\ ")
  



  
  # Current function is a truss, print if problem is continuous or discrete
  if isTruss(currentFunctionName):
      print("\\textbf{{{}}}: {} \\\\".format("Case", trussCase))


  for i in algorithmsMaxEval:
    print("\\textbf{{{}}}: {} \\\\".format(i, algorithmsMaxEval[i]))
    # print("\\\ ")
    
  print("\\end{table}")
  print()

  # Function analysed successefully
  analysisTable.currentFunction +=1

  # sys.exit("obrigado")

def renameColumnsTitles(columnsTitles):
  for key, title in enumerate(columnsTitles):
    if title == "min":
      columnsTitles[key] = "Best"
    if title == "median":
      columnsTitles[key] = "Median"
    if title == "mean":
      columnsTitles[key] = "Average"
    if title == "std":
      columnsTitles[key] = "St.Dev"
    if title == "max":
      columnsTitles[key] = "Worst"
    if title == "size":
      columnsTitles[key] = "fr"
  return columnsTitles

def renameProblem(problemName):
  if 'Problem110' in problemName:
    problemName = "10 bar truss"
  elif 'Problem125' in problemName:
    problemName = "25 bar truss"
  elif 'Problem160' in problemName:
    problemName = "60 bar truss"
  elif 'Problem172' in problemName:
    problemName = "72 bar truss"
  elif 'Problem1942' in problemName:
    problemName = "942 bar truss"
  elif "Problem21" in problemName:
    problemName = "Tension|Compression spring design"
  elif problemName == "Problem22":
    problemName = "Speed reducer design"
  elif problemName == "Problem23":
    problemName = "Welded beam design"
  elif problemName == "Problem24":
    problemName = "Pressure vesel design"
  elif problemName == "Problem25":
    problemName = "Cantilever beam design"
  return problemName

def printTexTable(tableInfo, analysisTable):
  # Define small adjusts before starting printing

  # print("Tableinfo: {}".format(tableInfo))
  # sys.exit("saindo (printTexTable)")

  # By default, highlights the smallest values
  highlightedArr = np.nanmin(tableInfo["np2dArr"], axis=0)
  if tableInfo["highlightOption"] == "Max":
    highlightedArr = np.nanmax(tableInfo["np2dArr"], axis=0)

    
  # sys.exit("obrigado")
  tabAlign = "{{@{{}} l | {} @{{}}}}".format(tableInfo["texTableAlign"]*len(tableInfo["columns"])) 

  # Begin of table structure .tex
  print("\\begin{table}[h]")
  print("\\centering")
  print("\\resizebox{0.8\\textwidth}{!}{")
  print("\\begin{minipage}{\\textwidth}")
  print("\\caption{{{}}}".format(tableInfo["caption"]))
  # print("\\vspace{{{}}}".format("0.5cm"))
  print("\\begin{{tabular}} {}".format(tabAlign))
  print("\\hline")
  print("&", end=" ")
  print(*tableInfo["columns"], sep=" & ", end=" \\\ ") # Printing columnTitles
  print("\n\\hline")

  # print("np2dArr: {}".format(np2dArr))

  for rowIndex, row in enumerate(tableInfo["np2dArr"]):
    # Printing row titlesdiscor
    print(tableInfo["rows"][rowIndex], end=" & ")
    for columnIndex, item in enumerate(row):
      printableItem = item
      if tableInfo["showInPercentage"]:
        printableItem = (item / analysisTable.problemsSize) * 100
        
      # Verifies if current item from matrix is the minimum, for highlightning
      if np.round(tableInfo["np2dArr"][rowIndex][columnIndex], tableInfo["precRound"]) == np.round(highlightedArr[columnIndex], tableInfo["precRound"]):
        if columnIndex == len(row) -1:
          print("\\textbf{{{:{prec}}}}".format(printableItem, prec=tableInfo["precFormat"]), end=" \\\ ")
        else:
          print("\\textbf{{{:{prec}}}}".format(printableItem, prec=tableInfo["precFormat"]), end=" & ")
      else:
        # Last item of row (fr) its printed without floating points and uses line break
        if columnIndex == len(row) -1:
          print("{:{prec}}".format(printableItem, prec=tableInfo["precFormat"]), end=" \\\ ")
        else:
          print("{:{prec}}".format(printableItem, prec=tableInfo["precFormat"]), end=" & ")
    print()

  if tableInfo["footer"]:
    for text in tableInfo["footer"]:
      print("\\{}".format(text))


  print("\\end{tabular}")
  print("\\end{minipage}}")

  print("\\end{table}")

def printAnalysisTable(np2dArr, columnsTitles, rowsTitles, analysisTable): 

  # List containing minimum value of each column
  # minOfEachColumn = np.nanmin(np2dArr, axis=0)
  # tableInfo = {
  #   "np2dArr": np2dArr,
  #   "columns": columnsTitles,
  #   "rows": rowsTitles,
  #   "precRound": 4, # Precision round parameter for highlightning minimum value from column
  #   "precFormat": "e", # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  #   "caption": "Best perrforming technique in each mechanical engineering problem",
  #   "texTableAlign": "r", # Tex align (could be l, c, r and so on)
  #   "highlightOption": "max", # Max | Min - Highlight the max or the minimum in each column
  #   "showInPercentage": True, # True | False - True if show data biased its percentage
  # }
  # printTexTable(tableInfo);
  minOfEachColumn = np.nanmax(np2dArr, axis=0)

  # print("analysisTable {}".format(analysisTable))
  # print("analysisTable on printAnalysisTable: {}".format(rowsTitles))
  # Define print format, table aligmnent and round rpecision
  precRound = 4 # Precision round parameter for highlightning minimum value from column
  precFormat = "e" # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  if (type(precRound) == int):
    precFormat = ".{}f".format(precRound) # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation

  caption= "Best perrforming technique in each mechanical engineering problem"
  texTableAlign = "r " # Align on right (could be l, c, r and so on)
  tabAlign = "{{@{{}} l | {} @{{}}}}".format(texTableAlign*len(columnsTitles)) 

  # Begin of table structure .tex
  print("\\begin{table}[h]")
  print("\\centering")
  print("\\resizebox{0.8\\textwidth}{!}{")
  print("\\begin{minipage}{\\textwidth}")
  print("\\caption{{{}}}".format(caption))
  # print("\\vspace{{{}}}".format("0.5cm"))
  print("\\begin{{tabular}} {}".format(tabAlign))
  print("\\hline")
  print("&", end=" ")
  print(*columnsTitles, sep=" & ", end=" \\\ ") # Printing columnTitles
  print("\n\\hline")

  # print("np2dArr: {}".format(np2dArr))

  for rowIndex, row in enumerate(np2dArr):
    # Printing row titlesdiscor
    print(rowsTitles[rowIndex], end=" & ")
    for columnIndex, item in enumerate(row):
      printableItem = (item / analysisTable.problemsSize) * 100
      # Verifies if current item from matrix is the minimum, for highlightning
      if np.round(np2dArr[rowIndex][columnIndex], precRound) == np.round(minOfEachColumn[columnIndex], precRound):
        if columnIndex == len(row) -1:
          print("\\textbf{{{:{prec}}}}\%".format(printableItem, prec=".2f"), end=" \\\ ")
          # print("\\textbf{{{:{prec}}}}".format(item, prec=".0f"), end=" \\\ ")
        else:
          print("\\textbf{{{:{prec}}}}\%".format(printableItem, prec=".2f"), end=" & ")
          # print("\\textbf{{{:{prec}}}}".format(item, prec=".0f"), end=" & ")
      else:
        # Last item of row (fr) its printed without floating points and uses line break
        if columnIndex == len(row) -1:
          print("{:{prec}}\%".format(printableItem, prec=".2f"), end=" \\\ ")
          # print("{:{prec}}".format(item, prec=".0f"), end=" \\\ ")
        else:
          print("{:{prec}}\%".format(printableItem, prec=".2f"), end=" & ")
          # print("{:{prec}}".format(item, prec=".0f"), end=" & ")
    print()

  print("\\end{tabular}")
  print("\\end{minipage}}")

  print("\\\ ")


  print("\\end{table}")


def get_all_values(nested_dictionary):
  for key, value in nested_dictionary.items():
    if type(value) is dict:
      get_all_values(value)
    else:
      print("{} : {}".format(key, value))

def computeWilcoxonTest():
  # ORDEM DOS MÉTODOS
  # CMAES Equal HOF APM, CMAES Equal HOF DEB, 
  # CMAES Equal LF APM, CMAES Equal LF DEB
  # CMAES Linear HOF APM, CMAES Linear HOF DEB, 
  # CMAES Linear LF APM, CMAES Linear LF DEB
  # CMAES Superlinear HOF APM, CMAES Superlinear HOF DEB, 
  # CMAES Superlinear LF APM, CMAES Superlinear LF DEB

  resultsDict = {
    '10 bar truss continuous':{
      # Problem110c_CMAES Equal Hof + APM - m280000
      'CMAES Equal HOF APM': [5060.853660338485, 5076.6692998518065, 5060.853660338472, 5076.6692998518, 5060.853660338584, 5060.853660338469, 5060.853660338574, 5060.853660338486, 5060.853660338481, 5060.853660338474, 5060.853660338805, 5060.853660338471, 5060.853660338638, 5071.4856616619, 5060.853660338625, 5060.853660338786, 5060.853660338475, 5076.669299851828, 5060.853660338515, 5076.669299851796, 5060.853660338562, 5060.853660346927, 5060.853660338532, 5060.853660338543, 5060.853660338495, 5060.853660338513, 5076.669299851823, 5060.853660338504, 5060.85366033848, 5060.853660338474],
      # Problem110c_CMAES Equal Hof + DEB - m280000
      'CMAES Equal HOF DEB' : [5060.85366033847, 5076.669299851813, 5076.669299851822, 5060.853660338464, 5060.853660338966, 5076.669299851797, 5060.8536603386265, 5076.6692998518665, 5060.853660338466, 5060.853660338633, 5060.853660338891, 5076.6692998518265, 5060.8536603384855, 5076.669299851795, 5060.853660338598, 5076.6692998518, 5060.853660338471, 5060.853660338469, 5060.853660338475, 5060.853660338754, 5060.853660338489, 5060.85366033859, 5076.669299851903, 5060.853660338535, 5060.8536603384755, 5060.853660338569, 5060.853660338473, 5060.853660338507, 5076.669299851809, 5060.853660338464],
      # Problem110c_CMAES Equal Last factible + APM - m280000
      'CMAES Equal LF APM': [5060.853660338488, 5076.669299851809, 5060.853660338473, 5076.669299851804, 5060.853660338584, 5060.853660338473, 5060.85366033858, 5060.853660338491, 5060.853660338482, 5060.853660338481, 5060.85366033882, 5060.853660338474, 5060.85366033864, 5076.669299851827, 5060.853660338632, 5060.8536603388075, 5060.853660338476, 5076.669299851828, 5060.853660338534, 5076.669299851799, 5060.853660338563, 5060.85366034703, 5060.853660338535, 5060.853660338545, 5060.853660338498, 5060.853660338529, 5076.669299851823, 5060.853660338504, 5060.85366033848, 5060.853660338475],
      # Problem110c_CMAES Equal Last factible + DEB - m280000
      'CMAES Equal LF DEB': [5060.853660338471, 5076.669299851814, 5076.669299851823, 5060.853660338467, 5060.853660338971, 5076.669299851801, 5060.85366033863, 5076.669299851867, 5060.853660338469, 5060.853660338645, 5060.853660338891, 5076.669299851831, 5060.853660338489, 5076.6692998517965, 5060.853660338604, 5076.669299851807, 5060.853660338474, 5060.853660338469, 5060.8536603384755, 5060.853660338756, 5060.853660338493, 5060.853660338598, 5076.669299851903, 5060.8536603385555, 5060.85366033848, 5060.853660338573, 5060.853660338478, 5060.853660338509, 5076.669299851811, 5060.8536603384655],
      # Problem110c_CMAES Linear Hof + APM - m280000
      'CMAES Linear HOF APM': [5060.853660338495, 5060.853660338511, 5076.669299851828, 5060.853660338469, 5060.853660338466, 5060.853660338521, 5060.85366033848, 5060.853660338653, 5060.853660338625, 5060.853660338498, 5060.853660338548, 5060.853660338479, 5060.8536603384755, 5060.853660338485, 5060.853660338465, 5060.853660338461, 5060.8536603384655, 5076.669299851834, 5060.853660338468, 5060.853660340867, 5060.853660338478, 5076.66929985179, 5060.853660338481, 5060.853660338508, 5060.853660338503, 5060.853660338509, 5060.853660338479, 5060.853660338571, 5060.853660338479, 5060.8536603385],
      # Problem110c_CMAES Linear Hof + DEB - m280000
      'CMAES Linear HOF DEB': [5060.853660338478, 5060.853660338545, 5076.669299851791, 5060.8536603384755, 5060.853660338496, 5060.853660338505, 5060.853660338475, 5060.853660338916, 5060.853660338462, 5076.669299851801, 5060.8536603384655, 5060.85366033848, 5060.853660338482, 5060.853660338532, 5060.853660338626, 5076.6692998518365, 5060.85366033846, 5060.853660338505, 5076.669299851794, 5060.853660338598, 5060.853660338478, 5060.853660338615, 5060.853660338636, 5060.853660338465, 5060.853660338511, 5060.853660338482, 5076.669299851791, 5060.853660338468, 5076.669299851901, 5060.853660338486],
      # Problem110c_CMAES Linear Last factible + APM - m280000
      'CMAES Linearl LF APM': [5060.853660338495, 5060.853660338513, 5076.669299851828, 5060.853660338472, 5060.853660338468, 5060.853660338526, 5060.853660338482, 5060.853660338653, 5060.853660338625, 5060.853660338498, 5060.853660338548, 5060.85366033848, 5060.853660338478, 5060.853660338487, 5060.8536603384655, 5060.853660338463, 5060.853660338466, 5076.669299851836, 5060.853660338471, 5060.853660340962, 5060.853660338478, 5076.669299851791, 5060.853660338483, 5060.853660338508, 5060.853660338503, 5060.853660338509, 5060.853660338482, 5060.853660338573, 5060.853660338486, 5060.853660338503],
      # Problem110c_CMAES Linear Last factible + DEB - m280000
      'CMAES Linearl LF DEB': [5060.853660338478, 5060.853660338551, 5076.669299851794, 5060.853660338477, 5060.853660338496, 5060.853660338505, 5060.85366033848, 5060.853660338916, 5060.853660338463, 5076.669299851803, 5060.853660338467, 5060.853660338482, 5060.853660338482, 5060.853660338535, 5060.853660338626, 5076.669299851837, 5060.85366033846, 5060.853660338506, 5076.669299851796, 5060.8536603386, 5060.853660338478, 5060.8536603386165, 5060.853660338645, 5060.853660338466, 5060.853660338511, 5060.853660338484, 5076.669299851794, 5060.853660338471, 5076.669299851901, 5060.853660338487],
      # Problem110c_CMAES Superlinear Hof + APM - m280000
      'CMAES Superlinear HOF APM': [5060.853660338475, 5060.85366033876, 5060.853660338464, 5076.669299852242, 5060.853660338479, 5060.853660338491, 5060.853660338471, 5060.853660338465, 5076.6692998518065, 5060.853660338477, 5060.85366033849, 5060.853660338514, 5060.853660338472, 5076.669299851794, 5060.853660338666, 5060.853660338496, 5060.853660338475, 5060.853660338475, 5060.853660338483, 5060.853660338465, 5060.853660338508, 5060.853660338508, 5060.853660338488, 5060.853660338495, 5060.853660338574, 5060.853660338469, 5076.669299852509, 5060.853660338565, 5060.853660338469, 5076.669299851848],
      # Problem110c_CMAES Superlinear Hof + DEB - m280000
      'CMAES Superlinear HOF DEB': [5076.669299851826, 5060.853660338483, 5060.853660338472, 5060.853660338733, 5060.853660338498, 5060.853660338466, 5060.853660338473, 5076.66929985179, 5060.853660338562, 5060.853660338484, 5060.853660338646, 5060.853660338577, 5060.853660338477, 5060.853660338477, 5060.853660338719, 5060.853660338662, 5060.853660338497, 5076.669299852134, 5060.853660338482, 5060.853660338618, 5060.853660338473, 5060.853660338489, 5076.669299851795, 5060.853660338683, 5060.853660338522, 5076.669299851797, 5060.853660338526, 5060.853660338495, 5060.853660338482, 5060.853660338494],
      # Problem110c_CMAES Superlinear Last factible + APM - m280000
      'CMAES Superlinear LF APM': [5060.8536603384755, 5060.8536603387665, 5060.853660338469, 5076.669299852245, 5060.8536603384855, 5060.853660338493, 5060.853660338476, 5060.853660338466, 5076.669299851822, 5060.853660338483, 5060.8536603385, 5060.853660338514, 5060.853660338484, 5076.6692998518, 5060.853660338682, 5060.853660338502, 5060.853660338478, 5060.853660338478, 5060.853660338484, 5060.853660338471, 5060.853660338508, 5060.853660338508, 5060.853660338494, 5060.853660338496, 5060.853660338577, 5060.853660338477, 5076.669299852518, 5060.853660338568, 5060.85366033847, 5076.669299851857],
      # Problem110c_CMAES Superlinear Last factible + DEB - m280000
      'CMAES Superlinear LF DEB': [5076.669299851836, 5060.853660338491, 5060.853660338477, 5060.853660338736, 5060.853660338502, 5060.853660338469, 5060.853660338476, 5076.66929985179, 5060.853660338568, 5060.853660338489, 5060.853660338646, 5060.853660338581, 5060.853660338484, 5060.853660338491, 5060.853660338737, 5060.8536603386665, 5060.853660338499, 5076.669299852138, 5060.8536603384855, 5060.853660338626, 5060.853660338477, 5060.853660338495, 5076.669299851799, 5060.853660338688, 5060.853660338526, 5076.669299851801, 5060.853660338532, 5060.853660338504, 5060.8536603384855, 5060.8536603385055],
      # Érica (Default APM)
      'PSO Érica Default APM': [5061.511041, 5061.538650, 5062.504955, 5077.356029, 5078.905268, 5077.018098, 5061.146115, 5077.493202, 5061.195271, 5061.231886, 5061.715322, 5061.451856, 5079.313613, 5060.947307, 5061.398212, 5079.171137, 5077.158168, 5076.884838, 5076.722037, 5078.351565, 5064.552685, 5074.377489, 5060.965706, 5073.523867, 5078.085830, 5087.805227, 5061.650347, 5061.469250, 5061.711690, 5077.352487]
    },
    '10 bar truss discrete': {
      # Problem110d_CMAES Equal Hof + APM - m90000
      'CMAES Equal HOF APM': [6458.439657708293, 5739.237346821028, 5527.836203993039, 5621.856151999931, 5586.521658677632, 5637.760489598, 6127.93576230672, 5527.836203993039, 6551.16207360213, 5762.490360125411, 5636.368645883426, 5534.741606684525, 6025.805814044739, 5599.435632579044, 6277.89480017916, 5580.518230193663, 5665.059035321548, 5620.428567638678, 6499.989839684168, 5528.921658677633, 5638.399970687289, 5860.33072356698, 5598.497476548706, 6039.781659187811, 5540.923372919617, 5951.661866139881, 5735.3241260543955, 6063.764334026822, 6080.673865833774, 6147.755580075759],
      # Problem110d_CMAES Equal Hof + DEB - m90000
      'CMAES Equal HOF DEB': [5841.419346259832, 5998.38254154009, 5594.344463754501, 5589.222567791731, 5627.4242042991455, 5577.834126309485, 5681.219970177112, 5586.738749614551, 5719.620827298105, 6335.1072683211905, 5577.623684368079, 5623.678411914447, 5608.121658677633, 5892.729658473561, 5646.770645832408, 5815.994853192622, 5537.5378924935585, 5614.223347178152, 7488.019216021976, 5632.902904981657, 5798.275528592831, 6490.534021813093, 5535.036203993039, 5655.200827808282, 5655.200827808282, 5532.1210872636375, 6691.736774896855, 6010.0377365142385, 6212.406333975805, 5591.318230193663],
      # Problem110d_CMAES Equal Last factible + APM - m90000
      'CMAES Equal LF APM': [6459.519657708293, 5759.432775509068, 5535.036203993039, 6923.0756840619715, 5625.036203993039, 5659.338229683485, 6134.134021813094, 5548.456723413927, 6594.0989565666205, 5774.755632579043, 5641.2065419482315, 5558.856151999931, 6056.878749104374, 5607.764983175389, 6277.89480017916, 5580.518230193663, 5670.222567791731, 5631.947061369119, 6503.949839684168, 5537.718801607659, 5859.6644637545005, 5935.980827298104, 5618.6812427327795, 6039.781659187811, 5541.941606684526, 6136.48628147252, 5844.129658473561, 6078.848152204002, 6090.2491385424955, 6190.785709548347],
      # Problem110d_CMAES Equal Last factible + DEB - m90000
      'CMAES Equal LF DEB': [5841.419346259832, 6002.5490350664595, 5594.344463754501, 7548.296358952002, 8413.050905287762, 5602.774645730374, 5681.219970177112, 5673.288619631784, 6579.393137930281, 6407.163891830326, 5594.538229683485, 5630.427112852048, 5622.772360074394, 5943.641346718992, 5661.978749614552, 5826.542645526302, 5574.954126309485, 5614.223347178152, 7509.151423688297, 5661.922567281554, 6284.694852682445, 6491.0457100585245, 5535.036203993039, 6831.81765826949, 5655.200827808282, 5532.1210872636375, 6832.136462938216, 6080.001268984422, 6244.5432169402975, 5665.938749614551],
      # Problem110d_CMAES Linear Hof + APM - m90000
      'CMAES Linear HOF APM': [5586.1210872636375, 5635.717373072671, 5534.3162039930385, 6353.419735953043, 5831.853035474602, 5695.384879699353, 5525.758411914447, 5684.91139845701, 5566.118230193662, 5691.154177792414, 5669.86438601993, 5532.121087263637, 6213.695840041291, 5807.910749308445, 5579.47667142082, 5603.0795284907945, 5579.758411914448, 6621.308333414609, 5527.836203993039, 5585.203684878256, 6280.670385866874, 5615.220827298105, 5513.317373072669, 5511.358411914447, 5568.402827757264, 5633.600827808283, 5888.1059700240585, 5528.0869053897995, 5683.884385509751, 5674.044360278464],
      # Problem110d_CMAES Linear Hof + DEB - m90000
      'CMAES Linear HOF DEB': [5805.750437349805, 5573.256151999932, 5541.941606684526, 5576.305710568702, 5490.737892493558, 5557.938749614551, 5562.518230193663, 5581.603684878257, 5595.941606684526, 5574.7966714208205, 5617.029918439093, 5583.640152408074, 5533.656151999931, 5560.223347178152, 5912.060047911683, 5564.802827757264, 5593.731502953401, 5531.035632579043, 5509.717373072671, 5668.223347178152, 5585.203684878256, 5578.223347178152, 5610.954126309485, 5551.985424861706, 5513.423347178152, 5625.938178200556, 5520.954126309485, 5602.836723924105, 5560.820255884109, 6258.739735953042],
      # Problem110d_CMAES Linear Last factible + APM - m90000
      'CMAES Linearl LF APM': [5613.941606684526, 5635.717373072671, 5540.819970177112, 6379.955580075758, 5834.4391913008685, 5695.893996581808, 5525.758411914447, 5684.91139845701, 5566.118230193662, 5696.100619325678, 5683.42706136912, 7771.836982359104, 6215.239995918575, 5808.990749308443, 5579.47667142082, 8861.063476395653, 5586.958411914447, 6621.308333414609, 5540.819970177112, 5585.203684878257, 6820.834177792413, 5615.220827298105, 5513.317373072669, 6316.417113107138, 5568.402827757264, 5633.600827808283, 5896.837710262596, 5528.0869053897995, 7986.20259251284, 6227.714333261557],
      # Problem110d_CMAES Linear Last factible + DEB - m90000
      'CMAES Linearl LF DEB': [5805.750437349805, 8663.275941476615, 5541.941606684526, 5579.099425014758, 5490.737892493558, 5557.938749614551, 5562.518230193663, 5587.2065419482315, 6409.247891728291, 6079.010879036122, 5617.029918439093, 6065.104384999574, 5533.656151999931, 6598.282879240193, 5912.060047911683, 5634.79997068729, 5656.224204299145, 6143.711995102292, 6931.693553364957, 5695.2065419482315, 5620.118230193662, 5586.1210872636375, 5610.954126309485, 5551.985424861706, 5522.039580994079, 5630.266957382905, 5520.954126309485, 5602.836723924105, 5835.546281982699, 6276.584411251215],
      # Problem110d_CMAES Superlinear Hof + APM - m90000
      'CMAES Superlinear HOF APM': [5532.121087263637, 5541.941606684525, 5655.606541948231, 5513.317373072669, 5534.741606684525, 5620.035581096114, 5626.047788762434, 5776.241346718994, 5578.058931335336, 5523.256723413927, 5541.941606684526, 5591.683684878256, 5592.601087263638, 5498.374645730374, 5615.158749104374, 5638.118230193663, 5525.758411914447, 5568.676671420821, 5513.423347178152, 5564.802827757264, 5567.141606684525, 5559.657294827922, 5562.637373072671, 5531.035632579043, 5785.601787895132, 5649.137892493558, 5599.505424861705, 5752.343451164365, 5710.722411812412, 5692.92965847356],
      # Problem110d_CMAES Superlinear Hof + DEB - m90000
      'CMAES Superlinear HOF DEB': [5507.758411914447, 5534.741606684525, 5768.687242579726, 5513.423347178152, 5523.256723413927, 6447.5376320178475, 6542.659839939257, 5541.941606684525, 5507.758411914447, 5938.734852682444, 5667.941606684526, 5605.105710568703, 5548.456723413927, 5530.337892493559, 5725.986074010272, 5867.533606888596, 5647.259009069906, 5635.102282084733, 5560.223347178153, 5527.836203993039, 5490.737892493558, 5490.737892493558, 5513.423347178152, 5504.158411914446, 5532.238411914447, 5507.758411914447, 5577.199970687289, 5549.203684878256, 5655.606541948231, 5579.758411914448],
      # Problem110d_CMAES Superlinear Last factible + APM - m90000
      'CMAES Superlinear LF APM': [6168.53506065487, 5723.07171041565, 5739.958411914447, 5513.317373072669, 5534.741606684525, 5629.231477211937, 5626.047788762434, 5790.385502596278, 5584.4851911478145, 5530.456723413928, 5611.938749614551, 5591.683684878257, 5592.601087263638, 5511.358411914447, 5982.194645220196, 5667.941606684526, 5525.758411914447, 5568.676671420821, 7999.368956821711, 5566.118230193662, 5568.1210872636375, 5613.623684368079, 5562.637373072671, 5939.3216586776325, 5785.601787895132, 5702.856151999932, 6757.703684368079, 5752.343451164366, 5710.722411812412, 5692.92965847356],
      # Problem110d_CMAES Superlinear Last factible + DEB - m90000
      'CMAES Superlinear LF DEB': [5530.337892493559, 9664.356982359102, 5768.687242579726, 6810.989838663812, 5523.256723413927, 6448.977632017848, 6545.845034658318, 5541.941606684525, 5513.317373072669, 5938.734852682444, 5740.2362039930385, 5605.105710568703, 6094.3778147589865, 5530.337892493559, 5732.35132072244, 6940.817892493558, 5647.259009069906, 5800.206333975804, 5568.9584119144465, 5834.704384999573, 6042.005502086099, 5490.737892493558, 5513.423347178152, 5541.840775304997, 5574.235632579043, 5507.758411914447, 6581.6308270430145, 5549.203684878257, 5656.118230193663, 7556.090125901342],
      # Érica (Default APM) #nRun 180k
      # 'PSO Érica Default APM': [5638.399971, 5718.186723, 5509.717373, 5509.717373, 6347.034801, 5509.717373, 5509.717373, 5509.717373, 5509.717373, 5676.754178, 5638.399971, 5509.717373, 5638.399971, 5509.717373, 5514.402828, 5631.525607, 5528.086905, 5638.399971, 5676.754178, 5509.717373, 5528.086905, 5514.402828, 5528.086905, 5514.402828, 5676.754178, 5509.717373, 5509.717373, 5562.518230, 5739.887243, 5509.717373]
      # Érica (Default APM) #nRun 90k
      'PSO Érica Default APM': [5638.399971, 5540.819970, 5626.047789, 5509.717373, 5509.717373, 5509.717373, 5509.717373, 5629.241010, 5634.799971, 5626.047789, 5514.402828, 5689.986074, 5626.047789, 5514.402828, 5676.754178, 5683.884386, 5536.965191, 5638.399971, 5629.241010, 5514.402828, 5514.402828, 5509.717373, 5656.224204, 5509.717373, 5514.402828, 5540.406542, 6463.898099, 5638.399971, 5509.717373, 5626.047789]
    },
    '25 bar truss continuous': {
      # Problem125c_CMAES Equal Hof + APM - m240000
      'CMAES Equal HOF APM': [484.0514235664732, 484.05142356647485, 484.05142356647275, 484.05142356647633, 484.0514235664729, 484.0514235664728, 484.0514235664725, 484.05142356647633, 484.05142356647286, 484.0514235664725, 484.05142356647275, 484.05142356647275, 484.0514235664731, 484.0514235664751, 484.0514235664727, 484.05142356647264, 484.0514235664727, 484.05142356647366, 484.0514235664725, 484.05142356647264, 484.0514235664738, 484.05142356647275, 484.05142356647264, 484.05142356647264, 484.05142356647343, 484.0514235664726, 484.05142356647417, 484.0514235664741, 484.05142356647343, 484.05142356647326],
      # Problem125c_CMAES Equal Hof + DEB - m240000
      'CMAES Equal HOF DEB': [484.0514235664724, 484.05142356647696, 484.0514235664729, 484.05142356647224, 484.05142356647286, 484.0514235664737, 484.05142356647303, 484.0514235664726, 484.051423566474, 484.05142356647275, 484.05142356647303, 484.0514235664725, 484.0514235664724, 484.05142356647565, 484.0514235664725, 484.05142356647264, 484.0514235664727, 484.0514235664723, 484.0514235664729, 484.0514235664726, 484.051423566474, 484.0514235664739, 484.0514235664726, 484.05142356647343, 484.05142356647264, 484.05142356647434, 484.05142356647326, 484.0514235665049, 484.05142356647247, 484.0514235664726],
      # Problem125c_CMAES Equal Last factible + APM - m240000
      'CMAES Equal LF APM': [484.05142356647343, 484.051423566475, 484.051423566473, 484.05142356647644, 484.05142356647355, 484.05142356647315, 484.0514235664726, 484.05142356647656, 484.0514235664732, 484.0514235664728, 484.05142356647303, 484.0514235664735, 484.05142356647315, 484.05142356647565, 484.0514235664728, 484.0514235664728, 484.05142356647286, 484.05142356647394, 484.0514235664729, 484.05142356647275, 484.05142356647417, 484.05142356647343, 484.0514235664727, 484.0514235664732, 484.05142356647343, 484.0514235664731, 484.05142356647417, 484.0514235664744, 484.0514235664735, 484.0514235664736],
      # Problem125c_CMAES Equal Last factible + DEB - m240000
      'CMAES Equal LF DEB': [484.05142356647264, 484.0514235664772, 484.05142356647326, 484.0514235664724, 484.05142356647303, 484.05142356647406, 484.0514235664732, 484.05142356647303, 484.05142356647417, 484.051423566473, 484.0514235664737, 484.0514235664727, 484.0514235664729, 484.05142356647605, 484.0514235664729, 484.0514235664731, 484.05142356647315, 484.0514235664733, 484.0514235664737, 484.05142356647275, 484.0514235664747, 484.0514235664741, 484.0514235664729, 484.05142356647343, 484.051423566473, 484.05142356647445, 484.0514235664735, 484.05142356650515, 484.0514235664728, 484.051423566473],
      # Problem125c_CMAES Linear Hof + APM - m240000
      'CMAES Linear HOF APM': [484.0514235664722, 484.05142356647224, 484.0514235664767, 484.0514235664723, 484.0514235664728, 484.05142356647366, 484.0514235664762, 484.051423566482, 484.0514235664756, 484.0514235664752, 484.0514235664727, 484.05142356647536, 484.0514235664727, 484.0514235664728, 484.05142356647264, 484.0514235664725, 484.0514235664965, 484.05142356648474, 484.05142356647747, 484.05142356647247, 484.05142356648537, 484.0514235664727, 484.05142356647235, 484.05142356647224, 484.0514235664751, 484.0514235664887, 484.0514235664775, 484.05142356647286, 484.0514235664763, 484.0514235664733],
      # Problem125c_CMAES Linear Hof + DEB - m240000
      'CMAES Linear HOF DEB': [484.05142356647355, 484.0514235664731, 484.05142356647946, 484.05142356647383, 484.0514235664734, 484.05142356647286, 484.05142356647445, 484.0514235664728, 484.0514235664737, 484.0514235664764, 484.051423566473, 484.05142356647264, 484.0514235664733, 484.05142356647264, 484.05142356647536, 484.051423566475, 484.05142356648497, 484.05142356647485, 484.0514235664724, 484.0514235664723, 484.051423566474, 484.05142356648645, 484.0514235664724, 484.0514235664742, 484.05142356647275, 484.05142356647275, 484.0514235664756, 484.0514235664726, 484.0514235664727, 484.05142356647247],
      # Problem125c_CMAES Linear Last factible + APM - m240000
      'CMAES Linearl LF APM': [484.0514235664729, 484.05142356647247, 484.0514235664767, 484.0514235664727, 484.0514235664735, 484.0514235664737, 484.0514235664762, 484.0514235664823, 484.05142356647576, 484.05142356647536, 484.051423566473, 484.0514235664758, 484.051423566473, 484.0514235664736, 484.0514235664735, 484.0514235664731, 484.0514235664967, 484.05142356648486, 484.0514235664777, 484.0514235664726, 484.0514235664857, 484.0514235664732, 484.0514235664725, 484.05142356647235, 484.05142356647946, 484.0514235664887, 484.0514235664775, 484.05142356647315, 484.05142356647656, 484.0514235664734],
      # Problem125c_CMAES Linear Last factible + DEB - m240000
      'CMAES Linearl LF DEB': [484.05142356647445, 484.0514235664733, 484.0514235664796, 484.05142356647417, 484.05142356647355, 484.0514235664732, 484.05142356647474, 484.051423566473, 484.0514235664737, 484.0514235664768, 484.05142356647343, 484.0514235664731, 484.0514235664736, 484.05142356647303, 484.05142356647576, 484.051423566475, 484.0514235664851, 484.0514235664749, 484.0514235664729, 484.05142356647264, 484.0514235664743, 484.0514235664866, 484.05142356647366, 484.0514235664746, 484.0514235664731, 484.0514235664729, 484.0514235664761, 484.0514235664728, 484.0514235664728, 484.0514235664727],
      # Problem125c_CMAES Superlinear Hof + APM - m240000
      'CMAES Superlinear HOF APM': [484.0514235664726, 484.0514235664748, 484.0514235664747, 484.05142356647275, 484.05142356648264, 484.0514235664729, 484.051423566474, 484.0514235664761, 484.05142356648184, 484.05142356647286, 484.0514235664744, 484.0514235664733, 484.05142356647383, 484.0514235664754, 484.0514235664751, 484.05142356647434, 484.051423566473, 484.05142356648247, 484.0514235664723, 484.0514235664735, 484.0514235664727, 484.0514235664797, 484.0514235664818, 484.0514235664754, 484.05142356647497, 484.05142356647355, 484.05142356647366, 484.0514235664781, 484.05142356647787, 484.05142356647667],
      # Problem125c_CMAES Superlinear Hof + DEB - m240000
      'CMAES Superlinear HOF DEB': [484.05142356647366, 484.0514235664736, 484.05142356647406, 484.05142356647343, 484.05142356648264, 484.0514235664811, 484.0514235664743, 484.0514235664729, 484.0514235664952, 484.05142356647536, 484.05142356647474, 484.05142356647343, 484.0514235664852, 484.05142356648, 484.0514235664739, 484.051423566473, 484.0514235664734, 484.051423566474, 484.0514235664748, 484.0514235664752, 484.05142356647247, 484.0514235664728, 484.0514235664743, 484.0514235664725, 484.0514235664751, 484.05142356647315, 484.05142356647343, 484.05142356647366, 484.05142356647264, 484.05142356647315],
      # Problem125c_CMAES Superlinear Last factible + APM - m240000
      'CMAES Superlinear LF APM': [484.05142356647355, 484.0514235664755, 484.05142356647514, 484.0514235664745, 484.0514235664836, 484.05142356647514, 484.0514235664743, 484.05142356647696, 484.05142356648184, 484.05142356647355, 484.0514235664746, 484.05142356647394, 484.0514235664754, 484.05142356647616, 484.05142356647855, 484.0514235664748, 484.0514235664735, 484.051423566483, 484.05142356647275, 484.0514235664738, 484.051423566473, 484.0514235664806, 484.05142356648213, 484.0514235664754, 484.05142356647536, 484.05142356647485, 484.0514235664745, 484.0514235664796, 484.05142356647843, 484.05142356647696],
      # Problem125c_CMAES Superlinear Last factible + DEB - m240000
      'CMAES Superlinear LF DEB': [484.05142356647406, 484.0514235664749, 484.051423566475, 484.0514235664741, 484.0514235664834, 484.05142356648184, 484.0514235664746, 484.0514235664729, 484.0514235664962, 484.0514235664759, 484.0514235664751, 484.0514235664738, 484.05142356648616, 484.05142356719733, 484.0514235664744, 484.0514235664738, 484.0514235664736, 484.05142356647497, 484.05142356647605, 484.05142356647553, 484.0514235664745, 484.05142356647343, 484.05142356647485, 484.05142356647303, 484.0514235664757, 484.05142356647417, 484.05142356647366, 484.05142356647394, 484.0514235664729, 484.0514235664742]
    },
    '25 bar truss discrete': {
      # Problem125d_CMAES Equal Hof + APM - m20000
      'CMAES Equal HOF APM': [517.4752403358339, 504.75508708137846, 509.50802450600236, 500.55615731865737, 502.74646382756794, 528.5306222000179, 500.0007754544732, 549.2295519627393, 499.6997052171944, 516.5773808103918, 500.55615731865737, 504.75508708137846, 503.06262798926196, 500.97524033583386, 528.6783056263857, 511.5803182350154, 510.33600406420203, 501.53062220001794, 521.9878230184006, 506.8687882343707, 510.38138847229436, 509.7038407450251, 501.0050870813785, 515.313230271112, 509.8442926752206, 516.9255507613295, 535.4316627732919, 512.2252403358339, 507.3346298619207, 510.2293758636646],
      # Problem125d_CMAES Equal Hof + DEB - m20000
      'CMAES Equal HOF DEB': [518.7529344981725, 521.8102928464882, 509.5263105731127, 526.3104689455628, 528.2238485781955, 546.4709287089285, 564.5846123065642, 521.5846298619207, 518.0475340648469, 532.0320028627925, 536.6441472539357, 503.78278124371724, 499.6997052171944, 500.4198584716496, 507.7309326255137, 499.2507754544732, 501.6655468447443, 505.5816924372968, 499.2507754544732, 501.2763105731127, 543.3941472539356, 503.3601649805601, 535.9986053159103, 538.605835798298, 507.6953935902891, 555.9418427889516, 524.5860040642021, 522.5306222000179, 505.054783116376, 537.9758270881607],
      # Problem125d_CMAES Equal Last factible + APM - m20000
      'CMAES Equal LF APM': [519.1733917176305, 504.75508708137846, 509.50802450600236, 505.7166170820232, 504.0219989462076, 528.5306222000179, 500.0007754544732, 549.2295519627393, 500.22524033583386, 516.5773808103918, 521.3316924372969, 505.280622200018, 503.06262798926196, 501.72524033583386, 528.6783056263857, 512.3303182350154, 510.33600406420214, 501.53062220001794, 525.185974400197, 513.0930770166567, 511.49568544664356, 510.2293758636645, 501.0050870813785, 571.5642864986771, 509.8442926752206, 516.9255507613295, 557.2849338269234, 512.2252403358339, 507.86016498056017, 510.2293758636646],
      # Problem125d_CMAES Equal Last factible + DEB - m20000
      'CMAES Equal LF DEB': [521.9943112443619, 521.8102928464882, 511.6579247015524, 527.3615391828416, 565.1367396586888, 547.2209287089285, 564.5846123065642, 525.1101649805603, 518.7975340648469, 532.7820028627925, 717.4399324122762, 505.28278124371724, 540.4709287089285, 504.588941488826, 510.6783056263856, 499.4453935902891, 502.1910819633838, 533.2621424784763, 499.2507754544732, 501.2763105731127, 543.3941472539356, 504.6357000991996, 536.7486053159103, 538.605835798298, 507.6953935902891, 556.6655292893878, 525.3360040642021, 569.5062283884336, 506.0029344981725, 538.7258270881607],
      # Problem125d_CMAES Linear Hof + APM - m20000
      'CMAES Linear HOF APM': [512.7209287089285, 513.2464638275679, 518.5149148988892, 499.44539359028903, 501.3943233530101, 509.66830170980035, 512.5021561172479, 502.8943233530102, 502.7763105731127, 500.22524033583386, 516.2878712515471, 500.82738081039173, 501.4709287089285, 499.77631057311265, 500.55615731865737, 506.3615391828415, 501.0219989462074, 520.2676873193024, 499.44539359028903, 503.7676873193023, 501.08169243729685, 499.4453935902891, 502.8346298619207, 529.6229237622015, 526.2222918163463, 501.3943233530101, 501.3943233530101, 534.7358289789124, 504.424170098555, 516.0219989462074],
      # Problem125d_CMAES Linear Hof + DEB - m20000
      'CMAES Linear HOF DEB': [526.1478565989909, 501.14001172610494, 503.10722755593645, 500.94539359028903, 502.8346298619207, 519.5007754544733, 518.3957040157848, 526.1029159290313, 507.4766209986084, 500.22524033583375, 520.598604302126, 511.1144766074654, 510.9367703364785, 499.44539359028903, 500.49646382756794, 499.44539359028903, 500.82738081039173, 509.9581523719304, 501.0219989462076, 501.0050870813785, 516.5490732713392, 516.7892519143218, 502.6997052171944, 509.7230877526278, 504.14001172610494, 502.1741700985549, 501.14001172610494, 500.49646382756794, 502.32738081039173, 499.2507754544732],
      # Problem125d_CMAES Linear Last factible + APM - m20000
      'CMAES Linearl LF APM': [513.4709287089285, 513.2464638275679, 518.5149148988892, 499.44539359028903, 501.3943233530101, 509.66830170980035, 513.0930770166567, 504.3346298619207, 502.7763105731127, 500.22524033583386, 516.2878712515471, 500.82738081039173, 501.4709287089285, 500.3018456917521, 500.55615731865737, 506.36153918284157, 501.3943233530101, 520.2676873193024, 518.8910512855941, 503.7676873193023, 502.5816924372968, 532.6230998612762, 503.5846298619207, 529.6229237622015, 526.7478269349857, 501.3943233530102, 501.3943233530101, 534.7358289789124, 504.94970521719443, 516.0219989462074],
      # Problem125d_CMAES Linear Last factible + DEB - m20000
      'CMAES Linearl LF DEB': [526.1478565989909, 501.14001172610494, 503.10722755593645, 501.3943233530101, 502.8346298619207, 519.5007754544733, 521.5675418980175, 527.6029159290313, 508.00215611724786, 500.22524033583375, 522.098604302126, 513.2932224379417, 510.9367703364785, 499.4453935902891, 501.0219989462076, 559.5626104339051, 502.32738081039173, 510.4836874905699, 501.0219989462076, 554.7653366470051, 516.5490732713392, 519.8187575565811, 502.86800985344615, 510.67123913442424, 506.5881631079014, 502.1741700985549, 501.14001172610494, 500.49646382756794, 502.32738081039173, 500.22524033583375],
      # Problem125d_CMAES Superlinear Hof + APM - m20000
      'CMAES Superlinear HOF APM': [501.83169243729685, 500.22524033583386, 501.08169243729685, 499.6997052171944, 500.55615731865737, 504.0050870813786, 505.83600406420214, 503.088941488826, 514.7464638275679, 501.14001172610494, 505.5634063701865, 499.6997052171944, 504.6486349799154, 501.0050870813785, 503.088941488826, 503.390011726105, 501.1400117261049, 510.3346298619207, 499.9709287089285, 510.99646382756794, 501.1997052171943, 521.8012323150384, 503.39001172610494, 509.70799382821286, 503.80909474328126, 502.44539359028903, 508.6997052171943, 501.0050870813785, 504.411235217839, 501.8935449720856],
      # Problem125d_CMAES Superlinear Hof + DEB - m20000
      'CMAES Superlinear HOF DEB': [505.59860430212586, 508.3104689455627, 500.4198584716496, 499.6997052171944, 499.2507754544732, 500.49646382756794, 506.3432531157312, 500.2252403358338, 504.74215220066264, 501.0050870813785, 510.14001172610483, 506.2719989462074, 511.2564677441531, 501.1997052171943, 500.52631057311265, 525.5787439178096, 515.7209287089285, 507.5176873193023, 501.3943233530102, 499.2507754544732, 505.3104689455627, 504.7835596246418, 511.42052493947614, 506.7507754544732, 516.1527705077461, 501.0219989462074, 504.22955196273915, 500.49646382756794, 499.4453935902891, 499.2507754544732],
      # Problem125d_CMAES Superlinear Last factible + APM - m20000
      'CMAES Superlinear LF APM': [503.10722755593633, 500.22524033583386, 501.08169243729685, 499.6997052171944, 500.55615731865737, 504.530622200018, 505.83600406420214, 503.088941488826, 515.4964638275679, 501.14001172610494, 506.3134063701865, 499.6997052171944, 507.45401684409956, 501.0050870813785, 503.088941488826, 503.8090947432813, 501.1400117261049, 511.0846298619207, 499.9709287089285, 510.99646382756794, 501.1997052171943, 521.8012323150384, 505.0846298619207, 509.70799382821286, 652.7238495075264, 502.44539359028903, 509.4497052171943, 540.4556431559157, 504.411235217839, 501.8935449720856],
      # Problem125d_CMAES Superlinear Last factible + DEB - m20000
      'CMAES Superlinear LF DEB': [505.59860430212586, 508.3104689455627, 500.4198584716496, 499.6997052171944, 499.2507754544732, 500.49646382756794, 506.3432531157312, 500.2252403358338, 504.74215220066264, 501.0050870813785, 510.14001172610483, 506.2719989462074, 615.002843478588, 501.1997052171943, 501.0518456917521, 531.8525812862389, 570.8798141550884, 507.5176873193023, 501.3943233530102, 499.44539359028903, 505.3104689455627, 504.7835596246418, 511.42052493947614, 508.2507754544732, 517.6755507613295, 501.0219989462074, 505.72955196273915, 506.9623054551179, 515.4497052171945, 499.2507754544732],
      # Érica (Default APM)
      'PSO Érica Default APM': [485.905250, 484.854179, 485.048797, 484.854179, 484.854179, 485.048797, 485.048797, 484.854179, 485.768951, 486.217880, 485.048797, 485.048797, 490.192345, 485.048797, 485.048797, 484.854179, 488.650938, 484.854179, 510.336915, 484.854179, 485.048797, 484.854179, 484.854179, 485.905250, 484.854179, 485.905250, 486.099868, 484.854179, 484.854179, 485.048797]
    },
    '60 bar truss continuous': {
      # Problem160c_CMAES Equal Hof + APM - m12000
      'CMAES Equal HOF APM': [312.3640075749378, 312.035738400269, 312.55202778678216, 310.36663399435497, 312.9635935083404, 311.36426962819087, 312.76239711166875, 309.11484603690707, 310.40670290512054, 324.83298926466534, 336.14913518910913, 319.8759606715534, 311.48365770365587, 312.8826355662967, 312.75780296390786, 311.28373399390534, 313.8594615420688, 310.56483108353444, 316.11568613426834, 311.6084320956725, 319.375057112268, 320.76934402520646, 310.5355584597289, 311.6303224465733, 318.3791777065958, 312.1451199268232, 312.7778699456206, 313.4357913266385, 311.6119056965526, 318.41883548622684],
      # Problem160c_CMAES Equal Hof + DEB - m12000
      'CMAES Equal HOF DEB': [314.1947350152066, 311.6984705507743, 320.6169384971868, 311.2985890838916, 310.0032156621056, 310.5253552862844, 311.18880486055446, 321.8586637945524, 315.91158921361955, 314.70174688688985, 310.26234661513746, 311.5171213453827, 340.70108378113645, 311.85730048522015, 317.40143469408054, 335.89722378450756, 312.278595930773, 321.8163075596375, 316.6638450356972, 312.12742430947213, 312.061365722701, 312.45991777649147, 315.0321052388282, 310.4483381131403, 312.9983352303563, 334.0227265264228, 322.1801249937825, 313.09538826490757, 313.4085433475438, 320.6958766246271],
      # Problem160c_CMAES Equal Last factible + APM - m12000
      'CMAES Equal LF APM': [312.42002135011614, 312.035738400269, 313.2354964048341, 310.54202075703876, 312.9635935083404, 311.40213993415244, 312.76239711166875, 309.17080511929555, 310.4090966622146, 324.83298926466534, 336.3965667558108, 319.9186346599183, 311.48365770365587, 312.89534026887594, 312.75780296390786, 311.3265497791363, 313.9085612966151, 310.57829813475854, 316.1890323514887, 311.6084320956725, 319.51596374072807, 320.76934402520646, 310.5955946334856, 311.6303224465733, 318.3791777065958, 312.1609532183476, 312.79472718854964, 313.47753432661165, 311.6408470019619, 319.2013660299965],
      # Problem160c_CMAES Equal Last factible + DEB - m12000
      'CMAES Equal LF DEB': [314.2344933737402, 311.72602616242983, 320.77759770060413, 311.4270013032205, 310.10159418914213, 310.5631115297766, 311.23039651029416, 321.8586637945524, 315.91158921361955, 314.70174688688985, 310.2765680685403, 311.5171213453827, 340.70108378113645, 311.94992027199606, 317.40143469408054, 335.9914415484452, 312.278595930773, 321.9951186425113, 316.6638450356972, 312.15871157353104, 312.076592661131, 312.6774780967341, 315.0321052388282, 310.4618271505812, 313.0389715762818, 334.0901207789924, 322.2435941056813, 313.9908199917708, 313.4707938699453, 320.6958766246271],
      # Problem160c_CMAES Linear Hof + APM - m12000
      'CMAES Linear HOF APM': [309.25770057549323, 309.3282587115802, 310.697957048511, 311.1768436671544, 317.95592455437804, 310.1660690794673, 309.77352428769365, 311.07409640592107, 313.84223605347717, 312.8777495137864, 310.5142937766801, 310.14820884492957, 309.6304097549322, 310.1685680509768, 310.37574701281767, 313.71988948117536, 313.2996003807903, 309.9070551671379, 311.0924479862485, 311.1786929030049, 310.1709136630266, 310.53683871499425, 310.88106031398956, 310.75658468898456, 309.93066390631174, 309.5376066948527, 309.61258903360215, 315.0769926107755, 310.11078111304107, 309.90188184940257],
      # Problem160c_CMAES Linear Hof + DEB - m12000
      'CMAES Linear HOF DEB': [312.23804209906854, 309.76581753312337, 309.42359897985074, 311.0158059854915, 311.14769289953483, 313.59571543255913, 310.64038768372467, 311.4969363070253, 310.4913644689282, 313.52172025619467, 311.4826046020287, 309.5120283999008, 310.0261531521557, 311.42665452033816, 309.6573879497809, 309.3603430676113, 309.5810616072622, 312.59453601098073, 309.4398624575488, 312.1673547582231, 314.7555465698025, 311.99056482413397, 311.5113723926857, 309.4686341586565, 310.49867175316353, 322.25903154532824, 313.33882905164177, 309.1902500901765, 310.4159779427818, 309.6988756434677],
      # Problem160c_CMAES Linear Last factible + APM - m12000
      'CMAES Linearl LF APM': [309.25918021781393, 309.3282587115802, 310.70683313636323, 311.1768436671544, 318.76832461206095, 310.1660690794673, 309.7745287186765, 311.07409640592107, 313.84392607180445, 312.8860621041125, 310.51532486321685, 310.14820884492957, 309.6397127764972, 310.17046594385516, 310.37620429088076, 313.73054557333, 313.3006019055473, 309.9070551671379, 311.0924479862485, 311.18909692791317, 310.19873061741254, 310.5564064768197, 310.89966259713145, 310.8197674408968, 309.94242310909414, 309.5376066948527, 309.61417766036953, 315.0983492673349, 310.11078111304107, 309.90188184940257],
      # Problem160c_CMAES Linear Last factible + DEB - m12000
      'CMAES Linearl LF DEB': [312.2394339245078, 309.7731773106804, 309.4247065875292, 311.0204155838427, 311.175830527148, 313.59571543255913, 310.6633178692759, 311.5221625521683, 310.4913644689282, 313.54749939867634, 311.508629475203, 309.57036415784137, 310.0261531521557, 311.4291784142761, 309.6618277606582, 309.3603430676113, 309.6066087424751, 312.59820239066306, 309.44683160093575, 312.16930232038163, 314.7744326510926, 311.99056482413397, 311.5692940995975, 309.47730731403686, 310.5216375581193, 322.3056278139477, 313.33882905164177, 309.19055883753794, 310.4159779427818, 309.70229507575624],
      # Problem160c_CMAES Superlinear Hof + APM - m12000
      'CMAES Superlinear HOF APM': [310.0347447733507, 310.1618260367493, 309.799749602734, 310.02642022438073, 309.3520163437224, 309.18622692428124, 309.71044173217905, 309.99031690639924, 309.5220421512347, 309.43358460191723, 311.4818447847424, 310.116869553311, 310.0221492094698, 309.92075615106376, 309.3480862326737, 309.6381815394379, 310.18326523622915, 311.6835053701992, 311.16921743610396, 309.19563311997564, 309.7386847013366, 311.8665834691739, 309.855654459651, 310.1695737097051, 310.12627012671896, 310.6906682045144, 309.75301137312545, 310.549932052282, 309.5025988477955, 309.8515376369496],
      # Problem160c_CMAES Superlinear Hof + DEB - m12000
      'CMAES Superlinear HOF DEB': [309.3868092355271, 309.6847832224485, 310.43716688764187, 310.0692276168227, 309.45118023459577, 310.89952308804124, 309.9120722072819, 310.53072517318697, 315.46292772318014, 310.27191611209815, 314.4258251069817, 309.920172382651, 309.3944375989514, 318.3470363570511, 309.35709806820614, 309.54666230933606, 310.35909938632136, 309.8063259548766, 309.69050176795776, 310.4671032507812, 312.06394492534963, 311.79962164631144, 311.37992801552355, 311.47092817838706, 310.03703609539855, 309.7343175799509, 310.82252538806034, 310.52530312087754, 309.4324719686644, 309.1212133190767],
      # Problem160c_CMAES Superlinear Last factible + APM - m12000
      'CMAES Superlinear LF APM': [310.03653061963416, 310.1645795966629, 309.80738676666954, 310.0276503068836, 309.3529374852751, 309.18622692428124, 309.7486644706838, 310.00043651451267, 309.5271603845729, 309.4353500265414, 311.4818447847424, 310.13220455572997, 310.06402128399117, 309.92075615106376, 309.3480862326737, 309.6384178973388, 310.18326523622915, 311.69078708334393, 311.17814697673117, 309.1960124348209, 309.7386847013366, 311.86715006492807, 309.85658504289364, 310.1737482062695, 310.5824603440159, 310.71180488953723, 309.7530784989439, 310.55592450129126, 309.5143931780037, 309.8533736504135],
      # Problem160c_CMAES Superlinear Last factible + DEB - m12000
      'CMAES Superlinear LF DEB': [309.38693716397665, 309.6847832224485, 310.44257270978017, 310.0810920099537, 309.4569937792459, 310.89952308804124, 309.94894933056474, 310.53072517318697, 315.46292772318014, 310.2837800866502, 314.4258251069817, 309.92034897432757, 309.4013150384468, 318.3470363570511, 309.35735102042463, 309.5634199177015, 310.37714606247454, 309.8063259548766, 309.6905041752316, 310.5169491563137, 312.0658998717208, 311.80057046921445, 311.43215770851583, 311.4828490674597, 310.03703609539855, 309.7343175799509, 310.952354215978, 310.5359680027868, 309.4698429717214, 309.1215743969607],
      # Érica (Default APM)
      'PSO Érica Default APM': [333.780272, 317.159422, 319.153539, 314.155900, 334.765591, 336.493622, 380.653120, 328.758153, 320.723661, 349.548837, 313.435021, 345.515483, 315.975486, 353.227225, 333.958327, 314.781617, 314.302051, 313.786156, 322.064294, 311.717441, 314.394537, 325.425763, 313.086119, 332.075586, 347.818036, 315.626156, 318.117280, 315.078139, 368.746395, 312.979341]
    },
    '72 bar truss continuos': {
      # Problem172c_CMAES Equal Hof + APM - m35000
      'CMAES Equal HOF APM': [379.6148111653862, 379.6148356208523, 379.61645514367706, 379.61523794625225, 379.6148135466999, 379.61482924231325, 379.6148032913488, 379.61481979181485, 379.61481183527843, 379.6148128802449, 379.6148496476542, 379.6148286786888, 379.61493615514195, 379.614817657441, 379.61481172467336, 379.6148162926263, 379.61494097253876, 379.6148431445661, 379.61485986345866, 379.61484349079575, 379.61488367222006, 379.6148104023905, 379.61575471884265, 379.6148593147044, 379.6148421531341, 379.6148071554202, 379.61482303429125, 379.6148477634671, 379.61480435501556, 379.61482425237097],
      # Problem172c_CMAES Equal Hof + DEB - m35000
      'CMAES Equal HOF DEB': [379.6148042758564, 379.61486254089635, 379.6148103594225, 379.6148077117703, 379.6148033576005, 379.61480445330926, 379.61481284975616, 379.6148573057572, 379.614863995308, 379.6148062380756, 379.61480496654025, 379.61480663191674, 379.6148038126308, 379.61480684048706, 379.6148107005167, 379.6148313530435, 379.614980209038, 379.6148233837901, 379.6148055425911, 379.6148639771671, 379.6148183023271, 379.61498598511463, 379.6148140604771, 379.614814922846, 379.6148158648516, 379.6149000482518, 379.61481163412253, 379.61482064042724, 379.6148426036169, 379.6148407404003],
      # Problem172c_CMAES Equal Last factible + APM - m35000
      'CMAES Equal LF APM': [379.61481136024213, 379.61484138712683, 379.6164647573934, 379.61523794625225, 379.6148136106506, 379.6148292739454, 379.6148033344376, 379.61481979181485, 379.6148126923342, 379.6148131455859, 379.61485059213436, 379.61483031461324, 379.6149430059777, 379.61481778335684, 379.614812109178, 379.61481641114085, 379.61494316577654, 379.6148435482212, 379.6148660820105, 379.61484629492026, 379.61488444571194, 379.6148104023905, 379.61575471884265, 379.6148593147044, 379.6148421531341, 379.6148084839588, 379.6148236595818, 379.6148477634671, 379.614804536618, 379.61482452118537],
      # Problem172c_CMAES Equal Last factible + DEB - m35000
      'CMAES Equal LF DEB': [379.6148043125683, 379.61486345978966, 379.6148105302559, 379.61480800373505, 379.6148034016797, 379.61480454440175, 379.61481549543026, 379.6148618128493, 379.614863995308, 379.6148063266029, 379.61480522526955, 379.61480663191674, 379.61480383814535, 379.61480696309576, 379.614810709923, 379.61483225127915, 379.614980209038, 379.6148233837901, 379.6148055425911, 379.6148639771671, 379.61481832058246, 379.61498598511463, 379.6148140604771, 379.6148150026734, 379.6148175383205, 379.6149427054614, 379.6148119478199, 379.61482072168377, 379.61484281718316, 379.6148407404003],
      # Problem172c_CMAES Linear Hof + APM - m35000
      'CMAES Linear HOF APM': [379.61480386254874, 379.61480262918616, 379.6148025360093, 379.61480293005155, 379.6148024801673, 379.61480246846037, 379.61480279214993, 379.61480255977006, 379.6148027831073, 379.61480250801173, 379.61480306889587, 379.6148026720494, 379.6148024751717, 379.6148047295022, 379.61480261902585, 379.614802750785, 379.61480276323806, 379.6148028572075, 379.61480671942206, 379.6148053550349, 379.61480255859624, 379.61480317020823, 379.6148236198349, 379.6148027495756, 379.6148047047372, 379.6148908330489, 379.6148026792563, 379.6148025348315, 379.61480287511307, 379.61480271586436],
      # Problem172c_CMAES Linear Hof + DEB - m35000
      'CMAES Linear HOF DEB': [379.61480266675875, 379.6148030470151, 379.614802607864, 379.61480249785615, 379.6148024625748, 379.6148025338226, 379.6148059121466, 379.6148025083453, 379.61480253375754, 379.61480246057874, 379.61480377890086, 379.6148025275933, 379.61480251431476, 379.6148026030984, 379.6148026447865, 379.6148153342424, 379.6148028878041, 379.614802646264, 379.61480269572667, 379.61480247880326, 379.61480255575896, 379.6148026610389, 379.6148032576847, 379.6148039047617, 379.61480254701274, 379.6148027150694, 379.6148025710301, 379.61480308214146, 379.61480344705774, 379.6148026413241],
      # Problem172c_CMAES Linear Last factible + APM - m35000
      'CMAES Linearl LF APM': [379.61480389641434, 379.6148026345403, 379.6148025406127, 379.6148029435281, 379.61480248114583, 379.6148024694785, 379.61480279214993, 379.61480256243243, 379.6148027831073, 379.61480250801173, 379.6148030853974, 379.61480267262897, 379.6148024751717, 379.6148047295022, 379.614802629269, 379.6148027689912, 379.61480276323806, 379.61480290532046, 379.61480671942206, 379.6148053550349, 379.6148025589177, 379.61480317020823, 379.61482818500076, 379.61480275035615, 379.6148047047372, 379.614893419379, 379.61480268319644, 379.61480253594345, 379.61480287511307, 379.61480271586436],
      # Problem172c_CMAES Linear Last factible + DEB - m35000
      'CMAES Linearl LF DEB': [379.61480268406774, 379.61480319907156, 379.614802617216, 379.61480249962005, 379.61480246326875, 379.61480253632595, 379.6148059121466, 379.6148025145227, 379.61480253669265, 379.61480246205133, 379.6148037963634, 379.61480253122346, 379.61480251431476, 379.61480260635057, 379.6148026505878, 379.6148153342424, 379.614802893396, 379.6148026618734, 379.61480269572667, 379.6148024877415, 379.61480255575896, 379.6148026641802, 379.6148032576847, 379.6148040914788, 379.614802549942, 379.6148027150694, 379.6148025710301, 379.61480308214146, 379.61480344705774, 379.6148026413241],
      # Problem172c_CMAES Superlinear Hof + APM - m35000
      'CMAES Superlinear HOF APM': [379.6148062220421, 379.6148027003879, 379.6148049979546, 379.6148031544494, 379.6148026800658, 379.61480473422625, 379.61480267993284, 379.61481022824034, 379.6148104765394, 379.61480457730056, 379.614803368381, 379.61480499500783, 379.6148026918479, 379.61480523057844, 379.6148025855264, 379.6148485303004, 379.6148025897036, 379.61480648651457, 379.61480467226556, 379.6149757031476, 379.6148025358138, 379.61480440430375, 379.6148073803452, 379.61480572660355, 379.6148029792805, 379.61480259387, 379.6148026368238, 379.61480641595506, 379.6148074840359, 379.61480311337124],
      # Problem172c_CMAES Superlinear Hof + DEB - m35000
      'CMAES Superlinear HOF DEB': [379.6148028624787, 379.61480343385625, 379.61480343205574, 379.61480696944676, 379.61480318100456, 379.614807092964, 379.6148057069439, 379.61480268957223, 379.6148053923865, 379.6148028704033, 379.61480302003406, 379.61480378720347, 379.6148079914191, 379.614809840233, 379.6148249796865, 379.61480283375914, 379.6148034005628, 379.6148026686303, 379.61481717667806, 379.61481578750886, 379.6148032128583, 379.6148031120181, 379.6148068991762, 379.6148045967924, 379.61480312407673, 379.61480358435233, 379.61480470531393, 379.61480882104297, 379.6148068957212, 379.6148031541867],
      # Problem172c_CMAES Superlinear Last factible + APM - m35000
      'CMAES Superlinear LF APM': [379.6148062220421, 379.6148027059545, 379.61480508989393, 379.61480315967754, 379.61480268257793, 379.61480551267925, 379.61480270140066, 379.61481048623097, 379.61481047815556, 379.61480457990166, 379.614803492469, 379.61480500239657, 379.614802701208, 379.6148053537378, 379.61480259974076, 379.61485044817283, 379.61480259060755, 379.61480654204195, 379.6148046962919, 379.6149758307402, 379.6148025358138, 379.6148044048018, 379.6148073837055, 379.61480601291146, 379.6148029792805, 379.61480259937247, 379.61480264132695, 379.6148065139679, 379.614810969339, 379.61480311337124],
      # Problem172c_CMAES Superlinear Last factible + DEB - m35000
      'CMAES Superlinear LF DEB': [379.61480293378764, 379.61480343385625, 379.6148034796238, 379.6148070193907, 379.61480323562284, 379.614807092964, 379.61480572630444, 379.61480269385163, 379.61480543583866, 379.61480287103416, 379.6148030384147, 379.61480512530943, 379.6148079914191, 379.61480984642316, 379.6148269035686, 379.61480283380064, 379.6148034005628, 379.6148026765696, 379.61481717667806, 379.6148159933763, 379.61480336342424, 379.61480312910084, 379.6148069110304, 379.61480463033473, 379.61480312407673, 379.6148035994304, 379.61480470997566, 379.61480882104297, 379.61480697983853, 379.6148031597541],
      # Érica (Default APM)
      'PSO Érica Default APM': [379.728931, 379.705903, 379.776847, 379.734723, 379.733414, 379.755768, 379.761440, 379.690778, 379.665661, 379.739904, 379.678288, 379.712825, 379.800782, 379.780126, 379.688733, 379.652909, 379.794904, 380.204020, 379.743921, 379.690851, 379.755814, 379.816027, 379.665037, 379.711024, 379.667606, 379.698767, 379.668775, 379.710024, 379.881695, 379.758661]
    },
    'Tension|Compression spring design': {
      # Problem21_CMAES Equal Hof + APM - m36000
      'CMAES Equal HOF APM': [0.01301559573787023, 0.012669811977093133, 0.012683609278753545, 0.012747442251211417, 0.012669637204886884, 0.012682139996060196, 0.012740432130075419, 0.012668192839712987, 0.012721946741825052, 0.012931623646015262, 0.012678774073050377, 0.01267354714976463, 0.01266530964481332, 0.012864677687341945, 0.012666612675261197, 0.012780888430941311, 0.01266534151034554, 0.01266639883769224, 0.012676704478904023, 0.012667042381749085, 0.012665307661855318, 0.012672182395843218, 0.012673806110654036, 0.012679778516140943, 0.012676784661386043, 0.012692327114019126, 0.012668640792902558, 0.01303373659459071, 0.012668423811586052, 0.01268496765701428],
      # Problem21_CMAES Equal Hof + DEB - m36000
      'CMAES Equal HOF DEB': [0.014380398902097665, 0.012672622744532732, 0.01282969899237414, 0.01267084379491383, 0.012670542339734697, 0.012679250501745192, 0.012691861053229456, 0.012675704327194106, 0.012708133686420682, 0.012683658021521068, 0.012677256490935904, 0.012700054261428645, 0.012665234720073706, 0.01277924920361204, 0.012705643448935575, 0.012736805659113165, 0.012852686883288413, 0.012672432462606681, 0.01322399818500954, 0.01275568753924376, 0.012668100462768504, 0.01267837970785457, 0.01277373389887099, 0.012688029831381698, 0.012684463129068652, 0.012807258038866225, 0.012712255552216606, 0.012668125651040441, 0.012665300227165803, 0.012723474972202319],
      # Problem21_CMAES Equal Last factible + APM - m36000
      'CMAES Equal LF APM': [0.013015596185359433, 0.012669812082174797, 0.012683609278753545, 0.012747442334436995, 0.012669641840948357, 0.0126821431032032, 0.012740433063325373, 0.012668197049246342, 0.012721946741825052, 0.012931623646015262, 0.012678774073050377, 0.012673554876002706, 0.01266530964481332, 0.012864963460719686, 0.012666612675261197, 0.012780888430941311, 0.01266534151034554, 0.01266639883769224, 0.012676704478904023, 0.01266769584820651, 0.012665307664496061, 0.012672193132418985, 0.012673878274828555, 0.012679780073132448, 0.012676788699819124, 0.012692389080566581, 0.012668640792902558, 0.01303373659459071, 0.01266842929346186, 0.01268496765701428],
      # Problem21_CMAES Equal Last factible + DEB - m36000
      'CMAES Equal LF DEB': [0.014380683514547337, 0.012672622744532732, 0.01282969899237414, 0.012670844234007135, 0.012670542975855972, 0.012679252628589105, 0.012691861053229456, 0.012675735298597102, 0.012708133686420682, 0.012683658275363718, 0.012677645845134342, 0.012700060230904265, 0.012665234720540802, 0.01277924920361204, 0.012705662701731575, 0.012736836200139889, 0.012853275023159279, 0.012672433858315575, 0.01322399818500954, 0.01275568753924376, 0.012668100462768504, 0.012678392032771262, 0.012775014657575643, 0.012688029831381698, 0.01268447276828805, 0.012807258038866225, 0.012712255552216606, 0.012668125651040441, 0.012665300235649518, 0.012723532255305524],
      # Problem21_CMAES Linear Hof + APM - m36000
      'CMAES Linear HOF APM': [0.012665372254293135, 0.012665238224730927, 0.012665243006463149, 0.012665307990733114, 0.012665232818054678, 0.012665243347880096, 0.012665312680665506, 0.012665233856727706, 0.012665241420256686, 0.012665236080571237, 0.012665233438477466, 0.012665233659356448, 0.012665232939889417, 0.012665232791402007, 0.012665232806945677, 0.012733222381289024, 0.012665251670831884, 0.012665233030539822, 0.012665234226435777, 0.012665282298091812, 0.012665233133472927, 0.01266524429373358, 0.012665635410462466, 0.012665239852694306, 0.012665299550870174, 0.012665232805365523, 0.01266523902680491, 0.01266587577258936, 0.012665232795391454, 0.01266525149974614],
      # Problem21_CMAES Linear Hof + DEB - m36000
      'CMAES Linear HOF DEB': [0.01266531327630881, 0.012665244315245057, 0.012665233637095067, 0.012665268022683985, 0.01266523423323186, 0.012665236901730131, 0.012665233266173127, 0.01266523289426935, 0.012665370183004333, 0.012665242111104102, 0.012665232869141386, 0.012665235261776923, 0.012665519425561912, 0.01266523280171585, 0.012665237825346636, 0.012665252191383716, 0.012665242615892771, 0.012665232833000604, 0.012665720938726902, 0.012665238999717153, 0.012665232855139796, 0.012665236951834616, 0.012665247137979818, 0.01266526526345121, 0.012665233238231217, 0.012665300794693102, 0.012665232964270443, 0.012665881232588137, 0.012665233773515168, 0.012665252484502294],
      # Problem21_CMAES Linear Last factible + APM - m36000
      'CMAES Linearl LF APM': [0.012665372605677924, 0.012665238279177251, 0.012665243006463149, 0.01266530820730977, 0.012665232818054678, 0.012665243361340292, 0.01266531382936262, 0.01266523395309891, 0.012665241420256686, 0.012665236080571237, 0.012665233438477466, 0.012665233659356448, 0.012665232939927827, 0.012665232791402015, 0.012665232806945776, 0.01273322238128903, 0.012665251670831884, 0.012665233030539822, 0.012665234226435777, 0.01266528419595768, 0.01266523313361976, 0.01266524429373358, 0.012665637919643539, 0.012665239852694306, 0.012665299550870174, 0.012665232805365523, 0.012665239037457835, 0.01266587577258936, 0.012665232795392105, 0.01266525149974614],
      # Problem21_CMAES Linear Last factible + DEB - m36000
      'CMAES Linearl LF DEB': [0.012665315717423608, 0.012665244317346997, 0.012665233637095067, 0.012665268022683985, 0.01266523423426625, 0.012665236928092398, 0.01266523327352723, 0.012665232895572337, 0.012665370183004333, 0.012665242382579177, 0.012665232869145346, 0.012665235262426975, 0.012665521970632976, 0.01266523280171585, 0.012665237825346636, 0.012665252191383716, 0.012665242615892771, 0.012665232833018968, 0.01266572152750074, 0.012665239015170688, 0.012665232855139796, 0.012665236951834616, 0.012665247139148605, 0.012665265317522195, 0.012665233301959422, 0.012665300794693102, 0.012665232964270669, 0.012665881232588137, 0.012665233774363809, 0.012665252543545258],
      # Problem21_CMAES Superlinear Hof + APM - m36000
      'CMAES Superlinear HOF APM': [0.012665233384115243, 0.012665413895553842, 0.012667857537458712, 0.012665485407898015, 0.012665232788324054, 0.012665233005574102, 0.012665233506568651, 0.012665253077051943, 0.012665236266511025, 0.012665240032891984, 0.012665236234809177, 0.012665232870287034, 0.012665232957056752, 0.012665233296733221, 0.012665232798550998, 0.01266523360288356, 0.012665236099574324, 0.012665234082828358, 0.012665234448287044, 0.012665232806946494, 0.01266524212649626, 0.012665234105672604, 0.012665241107448797, 0.012665232799472799, 0.012665232808771515, 0.012665236168838474, 0.012665370580872989, 0.012665232790349378, 0.012665232885711152],
      # Problem21_CMAES Superlinear Hof + DEB - m36000
      'CMAES Superlinear HOF DEB': [0.012665309735426756, 0.012665232788364223, 0.012665234141631983, 0.012665267584907361, 0.012665232819478147, 0.012665233731684686, 0.012665240778611172, 0.012665232793764764, 0.01266523317850254, 0.012665232904190835, 0.01266523283436462, 0.012665235520399847, 0.012665233203348506, 0.012665234495743166, 0.012665232829932165, 0.012665232846827664, 0.012665232789365159, 0.01266523279709916, 0.012665236179708255, 0.01266524224745577, 0.012665232788739237, 0.012665232940779746, 0.012665260776402646, 0.012665340943988895, 0.012665232862659912, 0.012665233175508476, 0.01266523310264544, 0.012665232994903465, 0.0126652327893095, 0.012676591888866473],
      # Problem21_CMAES Superlinear Last factible + APM - m36000
      'CMAES Superlinear LF APM': [0.012665233384115264, 0.012665420001434812, 0.012667857537458712, 0.012665485941237798, 0.012665232788324054, 0.012665233005574102, 0.012665233506905415, 0.012665253674309632, 0.012665236269962004, 0.012665240120639066, 0.012665236234809177, 0.012665232870287034, 0.0126652329593063, 0.012665233319621135, 0.012665232798752422, 0.01266523360288356, 0.012665236099574324, 0.012665234082828358, 0.012665234448287058, 0.012665232806946494, 0.01266524212649626, 0.012665234109324259, 0.012665241341755241, 0.012665232799472799, 0.012665232808771553, 0.012665236168838474, 0.012665370580872989, 0.012665232790349394, 0.012665232886698819],
      # Problem21_CMAES Superlinear Last factible + DEB - m36000
      'CMAES Superlinear LF DEB': [0.012665309890585157, 0.012665232788368804, 0.012665234150761068, 0.012665267656583936, 0.012665232819697156, 0.012665233731684686, 0.012665240886630015, 0.01266523279377165, 0.01266523318326277, 0.01266523290551418, 0.01266523283436462, 0.012665235520399847, 0.012665233203348506, 0.012665234495743166, 0.012665232830425173, 0.012665232846827664, 0.012665232789365159, 0.01266523279709916, 0.012665236179708255, 0.012665242255281433, 0.012665232788739237, 0.012665232940779746, 0.012665260786090704, 0.012665343089817035, 0.012665232862659912, 0.012665233175508476, 0.01266523310264544, 0.012665232994903465, 0.01266523278931599, 0.012676591888866476],
      # Érica (APM MED 3)
      'PSO Érica APM MED 3': [0.012719, 0.014517, 0.012719, 0.012719, 0.017406, 0.012704, 0.013267, 0.012701, 0.012841, 0.015391, 0.013555, 0.014192, 0.016506, 0.013671, 0.012677, 0.012847, 0.016661, 0.012702, 0.012679, 0.012943, 0.012719, 0.012705, 0.013331, 0.013502, 0.012703, 0.015268, 0.012714, 0.012666, 0.013280, 0.013615]
    },
    'Speed reducer design': {
      # Problem22_CMAES Equal Hof + APM - m36000
      'CMAES Equal HOF APM': [2996.3481649685314, 3363.87337961336, 2996.34816507846, 2996.3481649685305, 2996.3481649685305, 2996.3481649685305, 4664.4561522884105, 2996.348164968543, 3363.87337961336, 2996.3481649685305, 2996.3481649685305, 2996.3481649685527, 2996.348164968531, 2996.348164968563, 2996.348164968531, 5171.456479442574, 2996.348164968534, 2996.3481649685305, 2996.3481650071553, 2996.348164968536, 2996.3481649685336, 2996.348164973249, 2996.3481649685596, 2996.34816496853, 2996.3481649685573, 2996.348164968535, 3363.87337961336, 3032.854983724501, 2996.3481649685664, 2996.3481649685386],
      # Problem22_CMAES Equal Hof + DEB - m36000
      'CMAES Equal HOF DEB': [2996.3481649685305, 2996.3481649685323, 4664.456153013493, 2996.4048238997634, 2996.348164968538, 2996.3481649685314, 2996.3481649685323, 2996.3481649685314, 2996.34816496853, 2996.3481649685305, 2996.3481649685305, 2996.348164968535, 2996.348694422683, 2996.3481649685305, 2996.348164968532, 2996.348164968531, 2996.34816496853, 2996.348164969672, 2996.348165330713, 2996.3481649685714, 2996.386078072856, 2996.3481649700207, 2996.3481649685295, 2996.348164968559, 4664.456152288408, 2996.348164994326, 2997.4392042685845, 4424.4418986028095, 2996.3481649685573, 2996.3481649685305],
      # Problem22_CMAES Equal Last factible + APM - m36000
      'CMAES Equal LF APM': [2996.348164968532, 4193.420870451677, 2996.3481654116495, 2996.348164968531, 2996.348164968531, 2996.3481649685314, 4664.456152288411, 2996.3481649685496, 4664.456152288617, 2996.3481649685305, 2996.3481649685305, 2996.3481649685527, 2996.348164968532, 2996.348164968568, 2996.3481649685314, 5171.456479442575, 2996.3481649685355, 2996.3481649685314, 2996.3481650081667, 2996.3481649685364, 2996.348164968534, 2996.348165097175, 2996.3481649685596, 2996.3481649685305, 2996.3481649685573, 2996.348164968535, 3758.369036474305, 3359.329352647862, 2996.3481649726036, 2996.3481649685386],
      # Problem22_CMAES Equal Last factible + DEB - m36000
      'CMAES Equal LF DEB': [2996.3481649685314, 2996.3481649685327, 4664.456153129879, 3173.327526889827, 2996.3481649685395, 2996.3481649685323, 2996.3481649685336, 2996.3481649685323, 2996.3481649685305, 2996.3481649685305, 2996.3481649685305, 2996.3481649685377, 3173.327526886874, 2996.348164968531, 2996.3481649685323, 2996.3481649685314, 2996.348164968531, 2996.348164970102, 2996.3555947081286, 2996.348164968572, 3011.454746392081, 2996.3481649800483, 2996.348164968531, 2996.34816496856, 4664.456152288409, 2996.348164994326, 3173.3275268976386, 4424.4418986028095, 2996.3481649685577, 2996.3481649685314],
      # Problem22_CMAES Linear Hof + APM - m36000
      'CMAES Linear HOF APM': [2996.3481649685373, 2996.34816496853, 2996.3481649685314, 2996.3481649685295, 2996.3481649685295, 2996.348164968539, 2996.3481649685295, 2996.3481649685295, 2996.3481649685295, 2996.34816496853, 2996.3481649685327, 2996.3481649685305, 4913.461569621622, 2996.3481649685295, 2996.3481649685295, 2996.3481649685314, 2996.3481649685295, 2996.3481649685314, 2996.3481649686887, 2996.3481649685295, 2996.348164968536, 2996.348164968541, 2996.3481649685295, 2996.3481649685295, 2996.3481649685505, 2996.3481649685295, 2996.348164968542, 2996.3481649685295, 2996.3481649685336, 2996.34816496853],
      # Problem22_CMAES Linear Hof + DEB - m36000
      'CMAES Linear HOF DEB': [2996.3481649685295, 2996.3481649686, 2996.3481649685323, 2996.3481649685305, 2996.3481649685295, 2996.3481649685295, 2996.3481649685314, 2996.3481649685314, 2996.3481649685327, 2996.3481649699725, 2996.3481649685327, 2996.348164968533, 2996.3481649685295, 2996.3481649685295, 2996.3481649685295, 5171.456479442573, 2996.3481649685305, 2996.34816496853, 2996.3481649685327, 2996.3481649685295, 2996.3481649685295, 2996.348164969535, 2996.3481649685305, 2996.3481649685314, 2996.3481649685295, 2996.34816496853, 2996.34816496853, 2996.3481649685295, 2996.34816496853, 2996.34816496853],
      # Problem22_CMAES Linear Last factible + APM - m36000
      'CMAES Linearl LF APM': [2996.3481649685377, 2996.34816496853, 2996.3481649685314, 2996.3481649685305, 2996.34816496853, 2996.3481649685687, 2996.348164968534, 2996.34816496853, 2996.34816496853, 2996.3481649685305, 2996.348164968533, 2996.3481649685305, 4913.461569621623, 2996.34816496853, 2996.3481649685295, 2996.3481649685323, 2996.3481649685305, 2996.3481649685314, 2996.3481649686896, 2996.3481649685305, 2996.348164968536, 2996.348164968541, 2996.34816496853, 2996.34816496853, 2996.3481649685514, 2996.34816496853, 2996.348164968542, 2996.3481649685305, 2996.348164968535, 2996.34816496853],
      # Problem22_CMAES Linear Last factible + DEB - m36000
      'CMAES Linearl LF DEB': [2996.34816496853, 2996.3481649686546, 2996.3481649685323, 2996.348164968531, 2996.34816496853, 2996.34816496853, 2996.348164968532, 2996.3481649685336, 2996.3481649685336, 2996.3481649746827, 2996.348164968533, 2996.348164968533, 2996.34816496853, 2996.34816496853, 2996.34816496853, 5171.456479442574, 2996.3481649685314, 2996.34816496853, 2996.348164968535, 2996.3481649685295, 2996.3481649685295, 2996.348164969535, 2996.3481649685305, 2996.3481649685314, 2996.3481649685295, 2996.34816496853, 2996.3481649685305, 2996.34816496853, 2996.3481649685305, 2996.34816496853],
      # Problem22_CMAES Superlinear Hof + APM - m36000
      'CMAES Superlinear HOF APM': [3359.329352647842, 2996.3481649685314, 2996.3481649685787, 2996.348164968536, 2996.3481649685828, 2996.3481649685514, 2996.3481649685305, 2996.348164968533, 2996.3481649685327, 2996.3481649686037, 2996.3481649685355, 2996.3481649685305, 2996.3481649685323, 2996.348164968531, 2996.348164968531, 2996.3481649685314, 2996.3481649685305, 2996.348164968532, 3363.87337961336, 2996.348164968532, 2996.348164968531, 2996.3481649685295, 2996.3481649685305, 2996.348164968536, 2996.3481649685305, 2996.348164968532, 3340.099314935592, 2996.348164968532, 2996.3481649685305, 2996.348164968536],
      # Problem22_CMAES Superlinear Hof + DEB - m36000
      'CMAES Superlinear HOF DEB': [2996.348164968733, 2996.348164968555, 2996.348164968531, 2996.3481649888267, 2996.348164968532, 2996.3481649685314, 2996.3481649685323, 2996.3481649685314, 2996.3481649685314, 2996.348164982946, 2996.3481649685305, 2996.348164968605, 2996.3481649685305, 2996.348164968532, 2996.348164968704, 2996.348164968557, 2996.348164968531, 2996.3481649685345, 2996.3481649685314, 2996.34816496853, 2996.3481649685327, 2996.3481649685336, 2996.3481649685305, 2996.348164968531, 2996.34816496854, 2996.348164968535, 4913.461569621625, 2996.348164968544, 2996.3481649685314, 2996.348164968531],
      # Problem22_CMAES Superlinear Last factible + APM - m36000
      'CMAES Superlinear LF APM': [3359.3293526478437, 2996.348164968532, 2996.348164968588, 2996.348164968536, 2996.348164968585, 2996.3481649685523, 2996.3481649685555, 2996.348164968534, 2999.5398594092458, 2996.348164968608, 2996.348164968538, 2996.3481649685314, 2996.348164968537, 2996.348164968539, 2996.3481649685327, 2996.348164968532, 2996.3481649685423, 2996.3481649685355, 3758.369036474281, 2996.3481649685323, 2996.3481649685373, 2996.3481649685327, 2996.3481649685327, 2996.348164968536, 2996.3481649685327, 2996.348164968534, 4913.461569621639, 2996.348164968537, 2996.348164968532, 2996.3481649685395],
      # Problem22_CMAES Superlinear Last factible + DEB - m36000
      'CMAES Superlinear LF DEB': [2996.3481649687346, 2996.348164968555, 2996.3481649685323, 2996.3481649888267, 2996.3481649685823, 2996.3481649685345, 2996.348164968533, 2996.3481649685314, 2996.3481649685423, 2996.34816498295, 2996.3481649685327, 2996.348164968605, 2996.348164968533, 2996.3481649685345, 2996.3481649764694, 2996.3481649685614, 2996.348164968533, 2996.348164968537, 2996.3481649685345, 2996.348164968531, 2996.3481649685473, 2996.3481649685364, 2996.3481649685327, 2996.34816496854, 2996.34816496854, 2996.348164968536, 4913.461569621637, 2996.3481649685496, 2996.3481649685314, 2996.348164968534],
      # Érica (APM MED 3)
      'PSO Érica APM MED 3': [3005.707178, 3007.458896, 2996.385379, 2996.372989, 2996.381070, 999999999.000000, 2996.376415, 2996.371152, 2996.378956, 2996.375446, 3007.454465, 2996.360906, 2996.366961, 2996.387044, 2996.408534, 2996.376506, 3005.701992, 2996.398116, 3016.775330, 3005.716417, 2996.361055, 3007.457123, 2996.376496, 3007.481575, 3007.460028, 2996.367674, 2996.362598, 2996.356314, 2996.378728, 3016.777031]
    },
    'Welded beam design': {
      # Problem23_CMAES Equal Hof + APM - m320000
      'CMAES Equal HOF APM': [2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917903, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.38113411689179, 2.4336261925593736, 2.3811341168917943, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917908, 2.3811341168917894, 2.3811341168917894, 2.381158937532337, 2.3811341168980107, 2.38113411689179, 2.3811341168917894],
      # Problem23_CMAES Equal Hof + DEB - m320000
      'CMAES Equal HOF DEB': [2.3811341168917894, 2.3811341168917908, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.381134116891849, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917908, 2.3811341254158753, 2.3811341168917894, 2.38113411689179, 2.38113411689179, 2.38113411689179, 2.381134116891789, 2.3811341168917894, 2.381134116892174, 2.38113411689179, 2.3811341168923263, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.3811341168917894],
      # Problem23_CMAES Equal Last factible + APM - m320000
      'CMAES Equal LF APM': [2.3811341168917917, 2.3811341168917903, 2.38113411689179, 2.3811341168917894, 2.381134116891791, 2.3811341168917908, 2.3811341168917903, 2.3811341168917894, 2.3811341168917903, 2.3811341168917908, 2.3811341168917903, 2.433626247773902, 2.381134116891795, 2.381134116891791, 2.381134116891792, 2.3811341168917908, 2.3811341168917908, 2.3811341168917903, 2.3811341168917903, 2.3811341168917908, 2.38113411689179, 2.3811341168917908, 2.3811341168917894, 2.381134116891792, 2.381134116891791, 2.38113411689179, 2.381158937532337, 2.381134116898011, 2.381134116891791, 2.3811341168917903],
      # Problem23_CMAES Equal Last factible + DEB - m320000
      'CMAES Equal LF DEB': [2.3811341168917908, 2.3811341168917908, 2.38113411689179, 2.3811341168917903, 2.3811341168917903, 2.3811341168917903, 2.3811341168918503, 2.381134116891791, 2.3811341168917903, 2.381134116891791, 2.381134116891791, 2.3811341168917908, 2.3811341168917894, 2.381134116891791, 2.381134125415876, 2.3811341168917908, 2.3811341168917908, 2.38113411689179, 2.3811341168917903, 2.3811341168917908, 2.3811341168917903, 2.381134116892175, 2.3811341168917908, 2.3811341168923272, 2.3811341168917903, 2.3811341168917908, 2.38113411689179, 2.3811341168917903, 2.3811341168917908, 2.38113411689179],
      # Problem23_CMAES Linear Hof + APM - m320000
      'CMAES Linear HOF APM': [2.3811341168917894, 2.381134116891794, 2.38113411689179, 2.3811341168917903, 2.3811341168917894, 2.38113411689179, 2.3811341171262255, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.381134116891791, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168919635, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894],
      # Problem23_CMAES Linear Hof + DEB - m320000
      'CMAES Linear HOF DEB': [2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.381134116891789, 2.3811341168917894, 2.38113411689182, 2.381134116891789, 2.3811341168917894, 2.3811341168917894, 2.381134116891789, 2.3811341168917894, 2.3811341168917894, 2.381134116891798, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894],
      # Problem23_CMAES Linear Last factible + APM - m320000
      'CMAES Linearl LF APM': [2.3811341168917894, 2.381134116891795, 2.38113411689179, 2.3811341168917908, 2.3811341168917894, 2.3811341168917908, 2.3811341171262264, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917917, 2.3811341168917903, 2.3811341168917903, 2.38113411689179, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917903, 2.38113411689179, 2.381134116891964, 2.3811341168917894, 2.38113411689179, 2.3811341168917908, 2.3811341168917903, 2.38113411689179, 2.38113411689179, 2.38113411689179, 2.3811341168917894, 2.38113411689179],
      # Problem23_CMAES Linear Last factible + DEB - m320000
      'CMAES Linearl LF DEB': [2.3811341168917903, 2.3811341168917903, 2.3811341168917894, 2.3811341168917903, 2.3811341168917903, 2.3811341168917894, 2.38113411689179, 2.3811341168917908, 2.3811341168917908, 2.38113411689179, 2.3811341168917903, 2.3811341168917903, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917903, 2.3811341168918205, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.38113411689179, 2.3811341168917894, 2.381134116891798, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917903],
      # Problem23_CMAES Superlinear Hof + APM - m320000
      'CMAES Superlinear HOF APM': [2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.381134116891789, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.381134116891789, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.3811341168917948, 2.3811341168917894, 2.381134116891789, 2.3811341168917894, 2.3811341168940285, 2.3811341168917894, 2.38113411689179],
      # Problem23_CMAES Superlinear Hof + DEB - m320000
      'CMAES Superlinear HOF DEB': [2.38113411689179, 2.3811341168917894, 2.3811341168917903, 2.3811341168917894, 2.381134116891789, 2.3811341168917903, 2.3811341168917908, 2.3811341168917894, 2.381134116891789, 2.38113411689179, 2.381134116893376, 2.381134116891789, 2.3811341168917894, 2.381134116891789, 2.3811341168917894, 2.3811341168917894, 2.381134116891789, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917894, 2.381134116891794, 2.38113411689179, 2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.381134116891789],
      # Problem23_CMAES Superlinear Last factible + APM - m320000
      'CMAES Superlinear LF APM': [2.38113411689179, 2.3811341168917894, 2.3811341168917894, 2.3811341168917934, 2.3811341168917943, 2.3811341168917894, 2.3811341168917934, 2.3811341168917903, 2.3811341168917934, 2.3811341168917903, 2.3811341168917903, 2.3811341168917917, 2.381134116891791, 2.381134116891791, 2.38113411689179, 2.3811341168917926, 2.3811341168917903, 2.381134116891791, 2.381134116891791, 2.381134116891791, 2.3811341168917943, 2.3811341168917894, 2.3811341168917957, 2.3811341168917948, 2.3811341168917894, 2.38113411689179, 2.3811341168917934, 2.381134116894029, 2.381134116891791, 2.3811341168917917],
      # Problem23_CMAES Superlinear Last factible + DEB - m320000
      'CMAES Superlinear LF DEB': [2.3811341168917908, 2.3811341168917903, 2.3811341168917988, 2.381134116891791, 2.381134116891791, 2.3811341168917908, 2.3811341168917934, 2.381134116891795, 2.381134116891793, 2.3811341168917903, 2.381134116893399, 2.3811341168917903, 2.3811341168917894, 2.3811341168917903, 2.381134116891791, 2.3811341168917988, 2.381134116891791, 2.38113411689179, 2.3811341168917908, 2.381134116891796, 2.3811341168917903, 2.3811341168917908, 2.381134116891791, 2.3811341168917943, 2.381134116891797, 2.3811341168917908, 2.381134116891801, 2.38113411689179, 2.38113411689179, 2.3811341168917908],
      # Érica (APM MED 3) #nRun: 32k
      # 'PSO Érica APM MED 3': [2.567967, 3.623356, 2.390618, 2.411939, 2.512645, 2.821102, 3.549347, 4.127299, 3.590989, 2.673580, 2.381238, 2.383933, 2.529964, 3.587004, 2.381162, 3.019611, 2.482301, 2.664934, 2.381230, 2.425029, 2.381267, 2.381446, 2.381337, 2.831506, 2.381259, 3.241451, 2.381274, 3.183043, 2.459556, 2.788458]
      # Érica (APM MED 3) #nRun: 320k
      'PSO Érica APM MED 3': [2.537770, 3.019198, 3.005156, 2.381176, 2.381290, 2.990616, 2.381228, 2.521823, 2.381170, 2.525694, 2.797996, 2.665724, 2.381168, 3.317920, 2.389651, 2.381177, 2.381146, 2.422141, 3.350641, 2.381325, 2.653051, 2.381156, 2.961060, 2.381157, 2.669972, 2.382357, 2.381146, 2.504115, 3.412289, 3.585389]
    },
    'Pressure vesel design': {
      # Problem24_CMAES Equal Hof + APM - m80000
      'CMAES Equal HOF APM': [7332.841507759518, 9128.05074727785, 10325.056744708787, 6820.4100801400955, 6059.714335048436, 6207.08629725945, 6820.4100801400955, 6820.4100801400955, 6445.530422884877, 6410.086759869195, 6771.596850735377, 7332.84150775952, 6974.932897866693, 6144.289934961047, 6319.651334245538, 6820.4100801400955, 6771.596850735377, 7332.84150775952, 6414.497468883205, 6771.601274603081, 7544.49251792507, 6120.393274166187, 6069.5612705668855, 6820.4100801400955, 6773.364041800776, 7332.841507759518, 7323.681609741965, 6090.526201685864, 7273.511094325949, 6153.232705085039],
      # Problem24_CMAES Equal Hof + DEB - m80000
      'CMAES Equal HOF DEB': [6771.596850735377, 6410.086759869195, 7273.511094325949, 6771.596850735377, 6410.086759869195, 7274.013024311955, 6820.4100801400955, 6820.4100801400955, 6466.201617789899, 6820.4100801400955, 6820.410080140096, 7273.511094325949, 6820.410080140096, 6771.596850735377, 7332.841507759518, 6771.596850735377, 6410.086759869195, 6820.4100801400955, 7888.37129065246, 6820.4100801400955, 6820.410080140096, 7273.511094325949, 6820.4100801400955, 6820.410080140096, 6771.596850735377, 6771.596850735377, 6820.4100801400955, 7273.511094325949, 7273.511094325949, 6820.4100801400955],
      # Problem24_CMAES Equal Last factible + APM - m80000
      'CMAES Equal LF APM': [7332.84150775952, 9128.05074727785, 12944.17610875595, 6820.410080140097, 6059.714335048437, 9792.831122142288, 6820.410080140097, 6820.410080140097, 18133.339544171038, 6410.086759869195, 6771.596850735377, 7332.841507759522, 8582.612056355942, 11766.734450814485, 7888.975298643222, 6820.410080140097, 6771.596850735377, 7332.841507759522, 6895.217685267709, 10637.719551260225, 7544.49251792507, 9068.655747672909, 6090.526201685867, 6820.4100801400955, 6820.410080140097, 7332.841507759524, 28010.97872177441, 6090.526201685864, 7273.511094325949, 6424.066289451608],
      # Problem24_CMAES Equal Last factible + DEB - m80000
      'CMAES Equal LF DEB': [6771.596850735377, 6410.0867598691975, 7273.51109432595, 6771.596850735377, 6410.086759869196, 7332.84150775952, 6820.410080140097, 6820.4100801400955, 6686.112719395651, 6820.410080140096, 6820.410080140098, 7273.511094325949, 6820.410080140097, 6771.596850735379, 7332.84150775952, 6771.596850735377, 6410.086759869195, 6820.410080140097, 8017.281636477693, 6820.4100801400955, 6820.410080140097, 7273.511094325949, 6820.4100801400955, 6820.410080140097, 6771.596850735379, 6771.596850735377, 6820.410080140097, 7273.511094325949, 7273.511094325949, 6820.410080140097],
      # Problem24_CMAES Linear Hof + APM - m80000
      'CMAES Linear HOF APM': [7369.537073921283, 6090.923915040256, 7273.511094325949, 6416.244238399841, 9927.386361714374, 8901.479151903786, 6166.744260853879, 7050.673234290167, 6772.002483865616, 7175.523163816381, 8092.903444752486, 7273.511094325949, 6410.086759869195, 6820.4100801400955, 6821.916634073346, 6090.526201685864, 7078.10418175123, 9208.491400922574, 6135.824779081345, 6820.4100801400955, 6089.892868776244, 6994.9145191880525, 6820.4100801400955, 6410.086759869195, 11960.910624989529, 6948.483546055386, 7466.120704401314, 6459.3963692442885, 6820.410080140096, 6410.086759869195],
      # Problem24_CMAES Linear Hof + DEB - m80000
      'CMAES Linear HOF DEB': [7332.841507759518, 6410.086759869195, 7544.49251792507, 6410.086759869195, 7544.49251792507, 6410.086759869195, 6820.4100801400955, 6820.4100801400955, 6771.596850735377, 6410.086759869195, 7273.511094325949, 6820.4100801400955, 6771.596850735377, 6771.596850735377, 6820.4100801400955, 6771.596850735377, 6370.7797127298545, 6771.596850735377, 6771.596850735377, 7332.841507759518, 6820.4100801400955, 7332.841507759518, 6410.086759869195, 6820.4100801400955, 7273.511094325949, 6452.410558066805, 6820.4100801400955, 6727.509790774897, 6771.596850735377, 6820.4100801400955],
      # Problem24_CMAES Linear Last factible + APM - m80000
      'CMAES Linearl LF APM': [7860.384059275722, 7040.12197405708, 7273.511094325949, 8867.0651363364, 13228.36920869924, 595130.5, 13776.953610396142, 7050.673234290167, 21441.969957654914, 8592.890545273965, 8726.598474844594, 7273.511094325949, 6410.086759869195, 6820.4100801400955, 6939.837389397559, 6090.526201685864, 7273.511094325949, 11153.567289032933, 9264.526228108152, 6820.4100801400955, 8017.28163647769, 8962.513384376743, 6820.410080140097, 6410.086759869195, 11960.910624989529, 10262.949341435218, 7582.218201764036, 6771.596850735377, 6820.410080140097, 6410.086759869195],
      # Problem24_CMAES Linear Last factible + DEB - m80000
      'CMAES Linearl LF DEB': [7332.84150775952, 6410.086759869195, 7544.49251792507, 6410.086759869195, 7544.49251792507, 6410.086759869195, 6820.4100801400955, 6820.4100801400955, 6771.596850735377, 6410.086759869195, 7273.511094325949, 6820.4100801400955, 6771.596850735377, 6771.596850735378, 6820.4100801400955, 6771.596850735378, 6370.779712729856, 6771.596850735378, 6771.596850735377, 7332.841507759518, 6820.4100801400955, 7332.84150775952, 6410.086759869195, 6820.4100801400955, 7273.511094325949, 6686.112719395651, 6820.4100801400955, 6771.596850735377, 6771.596850735377, 6820.4100801400955],
      # Problem24_CMAES Superlinear Hof + APM - m80000
      'CMAES Superlinear HOF APM': [7194.714317723079, 6736.821935248154, 7332.841507759518, 6410.086759869195, 6509.3500085222895, 6410.086759869195, 7273.511094325949, 7332.841507759518, 6410.086759869195, 7794.932689325999, 7050.673234290168, 7340.0484289661545, 7919.766047898006, 7050.673234290167, 8952.720603931826, 6771.596850735377, 6820.4100801400955, 7937.4611030413, 6880.649560587168, 6962.636035392701, 6771.596850735377, 7273.511094325949, 6771.596850735377, 7332.841507759518, 8896.420505337665, 6686.112719395649, 7973.122399721253, 6059.714335048436, 7551.830992481943, 7283.423136304851],
      # Problem24_CMAES Superlinear Hof + DEB - m80000
      'CMAES Superlinear HOF DEB': [6410.086759869195, 6410.086759869195, 7361.109016741759, 6413.315525800433, 6820.4100801400955, 6410.086759869195, 6771.596850735377, 6771.596850735377, 6771.596850735377, 6090.526201685864, 7332.841507759518, 7544.49251792507, 6771.596850735377, 6820.4100801400955, 6771.596850735377, 6820.4100801400955, 6410.086759869195, 6820.4100801400955, 6771.596850735377, 7056.198935504788, 6771.596850735377, 7332.841507759518, 6820.4100801400955, 6771.596850735377, 6771.596850735377, 7332.841507759518, 7332.84150775952, 6410.086759869195, 7332.841507759518, 6410.086759869195],
      # Problem24_CMAES Superlinear Last factible + APM - m80000
      'CMAES Superlinear LF APM': [16194.49523046748, 27019.65130502183, 7332.84150775952, 6410.086759869195, 10417.737170702596, 6410.086759869195, 7273.511094325949, 7332.84150775952, 6410.086759869196, 10908.76463249453, 7050.673234290195, 7340.0484289661545, 7919.766047898006, 7050.673234290167, 8953.848538998296, 6771.596850735381, 6820.410080140097, 11437.82949874887, 8250.42234386793, 14646.819127219223, 6771.596850735381, 7273.511094325949, 6771.596850735379, 7332.84150775952, 18404.341510280818, 6686.112719395649, 8980.307504277833, 6059.714335048437, 8738.300772287459, 461835.525],
      # Problem24_CMAES Superlinear Last factible + DEB - m80000
      'CMAES Superlinear LF DEB': [6410.086759869195, 6410.086759869195, 7544.49251792507, 6686.112719395649, 6820.410080140096, 6410.086759869195, 6771.596850735381, 6771.596850735381, 6771.596850735381, 6090.526201685872, 7332.84150775952, 7544.49251792507, 6771.596850735381, 6820.410080140097, 6771.596850735381, 6820.410080140097, 6410.086759869196, 6820.4100801400955, 6771.596850735379, 7273.51109432595, 6771.596850735379, 7332.841507759524, 6820.410080140097, 6771.596850735381, 6771.596850735381, 7332.84150775952, 7332.84150775952, 6410.086759869195, 7332.84150775952, 6410.086759869196],
      # Érica (APM MED 3)
      'PSO Érica APM MED 3': [6820.410168, 6090.526364, 6771.596964, 7332.841559, 6059.714360, 6059.714372, 6820.410185, 6059.714387, 6424.066392, 7332.841553, 6370.779797, 6370.779814, 6059.714374, 6820.410161, 7544.492518, 6059.714388, 7332.841604, 6370.779804, 6090.526322, 6090.526297, 6059.714439, 6370.779776, 6771.596982, 6090.526428, 6370.779797, 6090.526272, 6059.714396, 6059.714452, 6820.410170, 7544.492518]
    },
    'Cantilever beam design': {
      # Problem25_CMAES Equal Hof + APM - m35000
      'CMAES Equal HOF APM': [75086.45583129591, 82808.83136741146, 73802.1542227268, 69705.2968044993, 66851.56890426991, 75658.53489823153, 75610.23941252292, 65387.32088526501, 64578.20525738556, 69062.30273255467, 69072.76188271679, 72207.90011663288, 81749.44855396089, 69277.6867842763, 65853.8164896951, 69362.56476884174, 69918.0589539694, 64578.1940169594, 69520.67570988057, 70643.43079587154, 70850.80625360245, 70707.00737526498, 69277.11171906767, 72094.244457647, 70503.24303168262, 76547.37278477316, 70643.43079587154, 65909.50335023385, 71192.24915901518, 64578.19401693747],
      # Problem25_CMAES Equal Hof + DEB - m35000
      'CMAES Equal HOF DEB': [69512.5192153577, 70646.48658633248, 68693.43079587154, 71540.13506854196, 72390.77968264517, 68201.35606893402, 68293.43079587154, 70974.52160439445, 69093.43079587154, 72793.43079587154, 70954.17752787998, 72184.87599453979, 70131.74528024763, 70707.00737500872, 66101.72150598564, 69493.43079587154, 70707.00737500875, 68644.09422337834, 70707.0073766287, 64697.56733700789, 64965.02546279448, 71093.43079587154, 69436.09307608289, 68597.13015195058, 77259.56591794259, 64965.013232641104, 72793.43080233409, 69893.43079587154, 64904.0661402249, 64578.19886745923],
      # Problem25_CMAES Equal Last factible + APM - m35000
      'CMAES Equal LF APM': [76343.64646310417, 84315.85781708997, 74213.01683605254, 73714.41829484983, 66856.90101933644, 75721.52980719555, 81578.31886832265, 65387.32088526501, 64578.21138682246, 69512.51921535858, 70979.73260391383, 72218.18901024797, 87351.37037929747, 69278.09324562397, 74448.91453270374, 69436.13249133147, 70960.14818333324, 64578.19401696506, 69753.42792724399, 70643.43079587158, 73508.84632787568, 70707.00737543964, 69277.11172393852, 77938.01335598134, 70503.24303168265, 76547.37278477316, 70643.43079587158, 68293.43079587161, 71207.41404282732, 64578.194016937494],
      # Problem25_CMAES Equal Last factible + DEB - m35000
      'CMAES Equal LF DEB': [69512.51921535772, 70647.23243884686, 68693.43079587161, 73852.31986861434, 73070.9711804807, 68201.35606893405, 68293.43079587162, 72393.43079587155, 69093.43079587154, 72793.43079587158, 77601.82300858197, 72184.87599453979, 70791.28043916333, 70707.00737500875, 73742.20901394467, 69493.43079587155, 70707.00737500875, 68644.0942233838, 70707.00737679914, 76212.23370068523, 64965.025975227036, 71093.43079587157, 69436.09361752002, 68597.13015195061, 79115.16422818306, 64965.013232641104, 72793.43080334658, 69893.43079587158, 64965.05683792619, 64578.20038971278],
      # Problem25_CMAES Linear Hof + APM - m35000
      'CMAES Linear HOF APM': [80169.7321406446, 72165.74455779938, 69757.80539134299, 77202.84745330944, 65194.11307342158, 64965.01323263433, 64998.66942583244, 64578.19401693837, 74469.42718812288, 71131.2625431154, 67133.54765761117, 65701.91052543269, 68147.25542040367, 64614.75950511002, 64965.01323263433, 68227.75457116541, 76877.13506024721, 74567.20430633654, 75980.58691812537, 82596.367623664, 76171.66757752805, 66126.17058378899, 73906.67109000655, 65443.43079587153, 70655.74849528619, 71227.21253896283, 64578.20193522208, 64965.01323264808, 82444.7220611139, 67815.29927592054],
      # Problem25_CMAES Linear Hof + DEB - m35000
      'CMAES Linear HOF DEB': [70707.0073750087, 70707.0073750087, 68693.43079587154, 69277.11170919191, 68201.35606893401, 68201.35606893401, 83673.37230763822, 69436.0902028796, 64578.19401693746, 69493.43079587152, 72745.45713875946, 69898.67993741015, 70293.43079587152, 68643.43079587152, 73093.43079587152, 70707.0073750087, 70707.01134681723, 68644.09422334678, 68597.13015195055, 68201.35606893401, 68644.09422334678, 64578.19401693746, 64578.19401693747, 68201.35606893401, 68201.35607409145, 64578.19401693746, 69032.82419704601, 64578.19401693746, 66351.68570425407, 64965.01323263433],
      # Problem25_CMAES Linear Last factible + APM - m35000
      'CMAES Linearl LF APM': [80171.32156252001, 72194.67063369876, 69940.63218394983, 77203.64133153985, 65194.37803642851, 64965.01323263434, 64998.66942583244, 64578.194016938956, 76694.560044457, 71131.59077747927, 67197.34478920298, 66372.05437026061, 68147.65144364316, 64615.09809186528, 64965.01323263433, 68497.221430912, 79679.55171867741, 96056.28102538682, 75980.58691812537, 94407.5711585336, 76306.00750379279, 66126.17058378899, 86614.70269103385, 65443.43079587154, 71692.06302843554, 71591.38605604533, 64578.2125012086, 64965.01323265551, 82525.54147967666, 67815.29927592054],
      # Problem25_CMAES Linear Last factible + DEB - m35000
      'CMAES Linearl LF DEB': [70707.00737500872, 70707.00737500872, 68693.43079587154, 69277.11170919193, 68201.35606893404, 68201.35606893404, 83673.37230763822, 69436.09020287963, 64578.19401693747, 69493.43079587154, 72745.45713875949, 70767.85685398638, 70293.43079587154, 68643.43079587154, 73093.43079587154, 70707.00737500873, 75079.46399496683, 68644.0942233468, 68597.13015195057, 68201.35606893402, 68644.0942233468, 64578.19401693747, 64578.19401693747, 68201.35606893402, 68201.35607409145, 64578.194016937494, 69493.43079587154, 64578.19401693747, 75933.64931375324, 64965.01323263434],
      # Problem25_CMAES Superlinear Hof + APM - m35000
      'CMAES Superlinear HOF APM': [65432.34323671717, 66780.19116981394, 73069.71822062657, 68861.2632448917, 92056.9735353402, 80111.18787756051, 81871.72001140993, 83273.11670259116, 70303.46495133893, 68693.43079587154, 75897.27614162886, 78893.38094816467, 78212.5349973691, 64578.1942608569, 70308.7022047421, 71281.74155427258, 73367.29866338088, 76625.03500475665, 81352.0994005474, 76602.0811878404, 68693.43081018161, 77199.80998992405, 72445.99709933206, 64965.01323263434, 79512.89635063312, 64578.19401693747, 68748.66670162656, 64965.01323263434, 74286.19684893278, 75404.92597597245],
      # Problem25_CMAES Superlinear Hof + DEB - m35000
      'CMAES Superlinear HOF DEB': [70707.00737500873, 70894.1964401381, 69512.51921535766, 68293.43079587152, 68597.13015195055, 68644.09422334678, 69043.43079587154, 64578.19401693746, 71044.29621753652, 72831.00289691528, 70707.0073750087, 64578.194016937494, 72184.87599453978, 64578.19401693747, 70707.0073750087, 64578.19401693747, 70503.2430316826, 68597.13015195057, 68293.43438023549, 70707.00737500873, 65772.85747513158, 68146.98390275135, 69436.09020287961, 69093.43079587154, 70274.61585483333, 72184.87599453976, 69415.2821919359, 70503.2430317569, 64578.19401693746, 69893.43079587154],
      # Problem25_CMAES Superlinear Last factible + APM - m35000
      'CMAES Superlinear LF APM': [65432.34323671717, 68696.108788886, 74631.7198236368, 69257.50416751232, 92056.9735353402, 80134.1231859491, 83328.03048971918, 103687.0666584709, 70503.29834537677, 68693.43079587155, 76033.43050281932, 78954.36419098663, 78264.33657994917, 64578.194295030815, 70323.07458767938, 73451.37563799412, 74258.40462124125, 90866.09462996319, 85625.33046479967, 80788.26578552127, 68693.43081646365, 77199.80998992405, 74495.91566717648, 64965.01323263435, 97437.28350365645, 64578.19401693748, 103724.5439886474, 64965.01323263435, 82121.6553445851, 75888.69054169723],
      # Problem25_CMAES Superlinear Last factible + DEB - m35000
      'CMAES Superlinear LF DEB': [70707.0073750088, 71293.43079587158, 69512.51921535771, 68293.43079587155, 68597.1301519506, 68644.09422334684, 69043.43079587157, 64578.194016937494, 71893.43079587154, 74523.0502892654, 70707.00737500878, 64578.19401693758, 72184.87599453983, 64578.1940169375, 70707.00737500873, 64578.1940169375, 70503.24303168265, 68597.13015195062, 68293.43657056271, 70707.00737500876, 65773.55186665329, 68146.98390275135, 69436.09020287963, 69093.4307958716, 71493.43079589188, 72184.87599453979, 96252.77488442011, 70503.24303179639, 64578.19401693747, 69893.43079587158],
      # Érica (APM MED 3)
      'PSO Érica APM MED 3': [65443.487527, 72746.361708, 69493.487857, 68643.509184, 76243.501901, 69493.513235, 69493.542202, 64965.096708, 69243.481005, 75143.537536, 71043.498192, 65444.624805, 64578.225361, 65443.468828, 71043.472764, 66393.527588, 65443.455991, 68643.469905, 73843.461556, 69993.482480, 69993.483969, 68597.191907, 69993.509615, 68643.513721, 65443.489083, 64623.923998, 70113.150205, 72293.461985, 68293.515745, 65443.528484]
    }
  }
  print("Wilcoxon test")
  for problemName, algorithm in resultsDict.items(): # the basic way
    # print("Problem: {}".format(problemName)) 
    print("{}".format(problemName)) # prints the problem
    print('Base algorithm\tAlgorithm\tp-value')
    for algorithmName, values in algorithm.items():
      for name, value in algorithm.items(): 
        # print("{}: {} ".format(algorithmName, values))
        if name == 'CMAES Superlinear HOF DEB':
          # print("{} - {}: Ranksums({},{}): {}".format(name, algorithmName, value, values, ranksums(value,values)))
          # print("{}\t{}\t{}".format(name, algorithmName, ranksums(value,values)))
          print("{}\t{}\t{}".format(name, algorithmName, ranksums(value,values)[1]))
    
    print()

      # print("{}: {} ".format(algorithmName, values)) # print results from each algorithm


  # # # BEGIN TEST BLOCK # # #
  # myDict = {
  #   'problem1': {
  #     'alg1': [1,2,3],
  #     'alg2': [4,5,6]
  #   },
  #   'problem2': {
  #     'alg1': [10,20,30],
  #     'alg2': [20,30,40]
  #   }
  # }
  # for problemName, algorithm in myDict.items(): # the basic way
  #   # print("Problem: {}".format(problemName)) 
  #   print("{}".format(problemName)) # prints the problem
  #   for algorithmName, values in algorithm.items():
  #     for name, value in algorithm.items(): 
  #       # print("{}: {} ".format(algorithmName, values))
  #       if name == 'CMAES Superlinear LF APM' or name == 'alg1':
  #         # print("{} - {}: Ranksums({},{}): {}".format(name, algorithmName, value, values, ranksums(value,values)))
  #         print("Ranksum\t{}\t{}\t{}".format(name, algorithmName, ranksums(value,values)))

  #     # print("{}: {} ".format(algorithmName, values)) # print results from each algorithm
  # # # END TEST BLOCK # # #


if __name__ == '__main__':
  if RUN_WILCOXON_TEST:
    computeWilcoxonTest()
    sys.exit("Executed Wilcoxon test");
  columnsTitles = ["Best", "Median", "Average", "St.Dev", "Worst", "fr", "Run Count"]
  analysisTable = Analysis(columns=columnsTitles)

  solutions, functions, functionsMaxEval = readResults()
  analysisTable.functionsName = functions
  analysisTable.functionsMaxEval = functionsMaxEval



  makeAnalysisNew(solutions, functions, analysisTable)


  tableInfo = {
    "np2dArr": analysisTable.table,
    "columns": analysisTable.columns,
    "rows": analysisTable.rows,
    "precRound": 4, # Precision round parameter for highlightning minimum value from column
    "precFormat": ".0f", # Precision format for priting table. Can be '.xf' or 'e' for scientific notation
    "caption": "Best performing technique in each {} problem".format(PROBLEMS_TYPE),
    "texTableAlign": "r", # Tex align (could be l, c, r and so on)
    "highlightOption": "Max", # Max | Min - Highlight the max or the minimum in each column
    "showInPercentage": False, # True | False - True if show data biased its percentage
    "footer": [],
  }

  # printAnalysisTable(analysisTable.table, analysisTable.columns, analysisTable.rows, analysisTable)

  # printTexTable(tableInfo, analysisTable)

  # print("analysisTable on main: {}".format(analysisTable.table))
  # print("analysisTable on main: {}".format(analysisTable.rows))




# Tension/compression spring design
# Speed reducer design
# Welded beam design
# Pressure vesel design
# Cantilever beam design

# *---*---*---*---*---*---*---*---*---* Pandas Personal Helper ---*---*---*---*---*---*---*---*---*---*---*
# print(grouped.index.names) # Return indexes names, in this case, only "Algorithm" is returned
# print(grouped.index) # Each index is each algorithm "CMAES+APM", "CMAES+DEB" ...
# print(grouped.index.values) # Returns name of all algorithms

# print(grouped.columns.levels) # Return a list of list containing all column levels, from top to bottom.
# print(list(grouped.columns.levels[1].values)) # Return all names of columns given a index
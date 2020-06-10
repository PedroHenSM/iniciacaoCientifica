import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
import os
import sys
from pathlib import Path

PROBLEMS_TYPE = "Other Engineering" # Trusses | Other Engineering
SELECTED_INDIVIDUAL = "Hof" # Last | Hof
TRUSS_CASE = "Continuous" # Continuous | Discrete

def readResults(individualToPick, problemsType):
  if individualToPick == "Hof":
    individualToPick = "Hall of fame"
  elif individualToPick == "Last":
    individualToPick = "Last individual"
  else:
    sys.exit("Individual to pick not defined.")

  algorithms = ["DE", "CMAES"]
  if problemsType == "Trusses":
    functions = [10, 25, 60, 72, 942]
  elif problemsType == "Other Engineering":
    functions = [21, 22, 23, 24, 25] # Other engineering problems
    # functions = [21] # Other engineering problems

  constraintHandlingMethods = ["DEB", "APM"]
  seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  , 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
  # seeds = [1, 2]

  basePath = Path(__file__).parent
  solutions = []
  for f in (functions): # Functions
    dataOfFunction = []
    for a in (algorithms): # Algorithms
      for p in (constraintHandlingMethods): # Constraint handling methods
        for s in (seeds): # Seeds
          tempData = []
          # filePath = (basePath / "../results/t{}/{}_f{}_p{}_s{}.dat".format(f, a, f, p, s)).resolve() # Trusses
          filePath = (basePath / "../results/f{}/{}_f{}_p{}_s{}.dat".format(f-20, a, f-20, p, s)).resolve()
          file = open(filePath)
          countFactibleInd = 0
          while True:
            buffer = file.readline()
            if individualToPick in buffer: # Find individual for data analysis
              hasFactibleSolution = True
              updateBestIndividual = False
              buffer = file.readline() # Read one more
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
                tempData.append("Problem{}_{}+{}".format(f, a, p)) # TODO Verify is this works
                # tempData.append("{}_f{}_p{}".format(a, f, p)) # Old line, works
                dataOfFunction.append(tempData)
              break
    # Each list from solutions contains the data for each function
    solutions.append(dataOfFunction)

  # print("Solutions:")
  # print(*solutions)
  # sys.exit("Printing final solutions, thanks!")
  return solutions, functions

def getExtraResults(np2dArr, rowsTitles, index, case, function):
  # Trusses problems
  if function == 110: # 10 bar truss
    if case == "Continuous":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 5060.87, 5060.92, 5061.98, 3.93e+00, 5076.70, "-"] 
      duvde = ["DUVDE*", 5060.85, float("NaN"), 5067.18, 7.94e+00, 5076.66, "-"] 
      apm = ["APM*", 5069.08, float("NaN"), 5091.43, float("NaN"), 5117.39, "-"]

      # Append name of algorihtms
      rowsTitles.append(smde[0])
      rowsTitles.append(duvde[0])
      rowsTitles.append(apm[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    elif case == "Discrete":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 5490.74, 5490.74, 5495.99, 1.13e+01, 5529.30, 666]
      duvde = ["DUVDE*", 5562.35, float("NaN"), 5564.90, 0.6, 5565.04, 666]
      apm = ["APM*", 5490.74, float("NaN"), 5545.48, float("NaN"), 5567.84, 666]

      # Append name of algorihtms
      rowsTitles.append(smde[0])
      rowsTitles.append(duvde[0])
      rowsTitles.append(apm[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  elif function == 125: # 25 bar truss
    if case == "Continuous":
      # Define algorithms and literature results
      smde= ["SMDE k=2*", 484.06, 484.07, 484.07, 0.0107, 484.10]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0) # Appends values of given list (remove first index, which is the string)
    elif case == "Discrete":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 484.85, 485.05, 485.44, 0.693, 487.13]
      duvde = ["DUVDE*", 485.90, float("NaN"), 498.44, 7.66e+00, 507.77]
      apm = ["APM*", 485.85, float("NaN"), 485.97, float("NaN"), 490.74]

      # Append name of algorihtms
      rowsTitles.append(smde[0])
      rowsTitles.append(duvde[0])
      rowsTitles.append(apm[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  elif function == 160: # 60 bar truss
    if case == "Continuous":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 308.94, 309.42, 309.49, 0.464, 311.21]
      duvde = ["DUVDE", 309.44, float("NaN"), 311.54, 1.46e+00, 314.70]
      apm = ["APM", 311.87, float("NaN"), 333.01, float("NaN"), 384.19]

      # Append name of algorihtms
      rowsTitles.append(smde[0])
      rowsTitles.append(duvde[0])
      rowsTitles.append(apm[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    elif case == "Discrete":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 312.73, 314.20, 315.12, 3.98e+00, 335.88]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  elif function == 172: # 72 bar truss
    if case == "Continuous":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 379.62, 379.63, 379.65, 0.0341, 379.73]
      duvde = ["DUVDE", 379.66, float("NaN"), 380.42, 0.572, 381.37]
      apm = ["APM", 387.04, float("NaN"), 402.59, float("NaN"), 432.95]

      # Append name of algorihtms
      rowsTitles.append(smde[0])
      rowsTitles.append(duvde[0])
      rowsTitles.append(apm[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    elif case == "Discrete":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 385.54, 386.81, 386.91, 1.05e+00, 389.21]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  elif function == 1942: # 942 bar truss
    if case == "Continuous":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 149932.00, 171218.50, 174369.63, 1.72e+04, 230139.00]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    elif case == "Discrete":
      # Define algorithms and literature results
      smde = ["SMDE k=2*", 153010.00, 181899.00, 181127.23, 1.73e+04, 216514.00]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  # Other Engineering Problems
  elif function == 21: # Tension/compression spring
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 0.012666, 0.012892, 0.013131, 6.28e-4, 0.015318, 50] 
    apmSpor = ["APM SPOR", 0.012667602164, float('NaN'), 0.013748492439, float('NaN'), 0.017093902154, float('NaN')]
    apmMed3 = ["APM MED 3", 0.01266, 0.01312, 0.01389 , 9.1731e-03, 0.01777, 35]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == 22: # Speed reducer
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 2996.3483, 2996.3495, 2996.3501, 7.45e-3, 2996.3599, 50] 
    apmSpor = ["APM SPOR", 2996.34850933205, float('NaN'), 2996.35243640334, float('NaN'), 2996.36609677358, float('NaN')]
    apmMed3 = ["APM MED 3", 2996.3622 , 2996.3780, 2999.6083, 3.4911e+01, 3016.7808, 35]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == 23: # Welded beam
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 2.38335, 2.92121, 2.99298, 2.02e-1, 4.05600, 50] 
    apmSpor = ["APM SPOR", 2.38113481849464 , float('NaN'), 2.58228221674671 , float('NaN'), 3.20898593483156, float('NaN')]
    apmMed3 = ["APM MED 3", 2.38114 , 2.43315, 2.67102, 2.0656e+00, 3.46638, 35]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == 24: # Pressure Vesel
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 6059.855, 6426.710, 6545.126, 1.24E+2, 7388.160, 50] 
    apmSpor = ["APM SPOR", 6059.73045731256 , float('NaN'), 6581.18398763114 , float('NaN'), 7333.93495942434, float('NaN')]
    apmMed3 = ["APM MED 3", 6059.7143 , 6370.7797, 6427.6676, 2.6221e+03, 7544.4925, 35]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  elif function == 25:
    # Define algorithms and literature results
    proposedAISGA = ["Proposed AIS-GA", 64834.70, 74987.16, 76004.24, 6.93E+3, 102981.06, 50] 
    # best median avg std worst
    apmSpor = ["APM SPOR", 64599.980343 , float('NaN'), 66684.276327 , float('NaN'), 72876.210779, float('NaN')]
    apmMed3 = ["APM MED 3", 64578.271, 68294.702, 71817.816, 1.0431e+05, 173520.325, 35]

    # Append name of algorihtms
    rowsTitles.append(proposedAISGA[0])
    rowsTitles.append(apmSpor[0])
    rowsTitles.append(apmMed3[0])

    # Append reults on np2darr
    np2dArr = np.append(np2dArr, [proposedAISGA[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmSpor[1:]], axis=0)
    np2dArr = np.append(np2dArr, [apmMed3[1:]], axis=0)
  else:
    sys.exit("Function not defined.")
  return np2dArr, rowsTitles

# Generates statistical measures and generate .tex table
def makeAnalysis(solutions, functions):
  # Read solutions and store on dataframe
  case="discrete"
  for key, solution in enumerate(solutions):
    df = pd.DataFrame.from_records(solution)
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
    df.drop([col for col in df.columns if "DesignVariables" in col], axis=1, inplace=True)
    df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
    df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)

    # Group by algorithm and calcualtes statistical measures
    grouped = df.groupby(["Algorithm"])
    grouped = grouped.agg(["min", "median", "mean", "std", "max", "size"]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works
    # grouped = grouped.agg(["min", "median", "mean", "std", "max", np.size]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works
    # print("grouped: {}".format(grouped))
    # Generates 2d np array from pandas grouped df
    np2dArr = grouped.to_numpy()
    # Get row and columns titles
    rowsTitles = list(grouped.index.values) 
    columnsTitles = list(grouped.columns.levels[1])
    # Get extra results (appending other algorithms on 2d np arr)
    np2dArr, rowsTitles = getExtraResults(np2dArr, rowsTitles, key, case, functions[key])
    # sys.exit("thx")
    # Generates .text table
    np2dArrToLatex(np2dArr, columnsTitles, rowsTitles)

    # break

def np2dArrToLatex(np2dArr, columnsTitles, rowsTitles):
  # List containing minimum value of each column
  # minOfEachColumn = np2dArr.min(axis=0) # works
  minOfEachColumn = np.nanmin(np2dArr, axis=0)

  # Define print format, table aligmnent and round rpecision
  precRound = 4 # Precision round parameter for highlightning minimum value from column
  precFormat = "e" # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  if (type(precRound) == int):
    precFormat = ".{}f".format(precRound) # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  
  texTableAlign = "r " # Align on right (could be l, c, r and so on)
  tabAlign = "{{@{{}} l | {} @{{}}}}".format(texTableAlign*len(columnsTitles)) 

  # Gets caption biased on string before underscore
  caption = rowsTitles[0].split("_")[0]

  # Removes underscore from algorithm names
  for i in range(len(rowsTitles)):
    splitted = rowsTitles[i].split("_")
    if len(splitted) == 1:
      rowsTitles[i] = splitted[0]
    elif len(splitted) == 2:
      rowsTitles[i] = splitted[1]
    else:
      sys.exit("Unexpected splitted size. Please verify method")

  # Begin of table structure .tex
  print("\\begin{table}[h]")
  print("\\centering")
  print("\\caption{{{}}}".format(caption))
  # print("\\vspace{{{}}}".format("0.5cm"))
  print("\\begin{{tabular}} {}".format(tabAlign))
  print("\\hline")
  print("&", end=" ")
  print(*columnsTitles, sep=" & ", end=" \\\ ") # Printing columnTitles
  print("\n\\hline")

  for rowIndex, row in enumerate(np2dArr):
    # Printing row titles
    print(rowsTitles[rowIndex], end=" & ")
    for columnIndex, item in enumerate(row):
      # Verifies if current item from matrix is the minimum, for highlightning
      if np.round(np2dArr[rowIndex][columnIndex], precRound) == np.round(minOfEachColumn[columnIndex], precRound):
      # if np2dArr[rowIndex][columnIndex]== minOfEachColumn[columnIndex]:
        # Last item of row, use line break 
        if columnIndex == len(row) -1:
          print("\\textbf{{{:{prec}}}}".format(item, prec=precFormat), end=" \\\ ")
        else:
          print("\\textbf{{{:{prec}}}}".format(item, prec=precFormat), end=" & ")
      else:
        # Last item of row, use line break
        if columnIndex == len(row) -1:
          print("{:{prec}}".format(item, prec=precFormat), end=" \\\ ")
        else:
          print("{:{prec}}".format(item, prec=precFormat), end=" & ")
    print()

  print("\\end{tabular}")
  print("\\\ ")
  print("\\textbf{{{}}}: {}".format("Indíviduo Selecionado", SELECTED_INDIVIDUAL))
  print("\\end{table}")
  # Begin of table structure .tex
  print()

if __name__ == '__main__':
  solutions, functions = readResults(SELECTED_INDIVIDUAL, PROBLEMS_TYPE)
  makeAnalysis(solutions, functions)


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
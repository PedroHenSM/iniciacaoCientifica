import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
import os
import sys
from pathlib import Path

def readResults(individualToPick):
  if individualToPick == "hof":
    individualToPick = "Hall of fame"
  elif individualToPick == "last":
    individualToPick = "Last individual"
  else:
    sys.exit("Individual to pick not defined.")

  algorithms = ["DE", "CMAES"]
  functions = [10, 25, 60, 72, 942]
  constraintHandlingMethods = ["DEB", "APM"]
  seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  , 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
  # seeds = [1, 2, 3, 4, 5]

  basePath =  Path(__file__).parent

  solutions = []
  for f in (functions): # Functions
    dataOfFuction = []
    for a in (algorithms): # Algorithms
      for p in (constraintHandlingMethods): # Constraint handling methods
        for s in (seeds): # Seeds
          tempData = []
          filePath = (basePath / "../results/t{}/{}_f{}_p{}_s{}.dat".format(f, a, f, p, s)).resolve()
          file = open(filePath)
          while True:
            buffer = file.readline()
            if individualToPick in buffer: # Find individual for data analysis
              buffer = file.readline() # Read one more
              buffer = buffer.split(" ")
              if p == "DEB":
                if float(buffer[1]) != 0:
                  sys.exit("Individual infactible. Problem: {}: {}+{}+s{}".format(f, a, p, s))
              elif p == "APM":
                if float(buffer[0]) != float(buffer[1]):
                  sys.exit("Individual infactible. Problem: {}: {}+{}+s{}".format(f, a, p, s))
              for item in buffer:
                tempData.append(float(item))
              
            elif "CPU time used" in buffer:
              buffer = file.readline()
              tempData.append(float(buffer))
              tempData.append("Problem{}_{}+{}".format(f, a, p)) # TODO Verify is this works
              # tempData.append("{}_f{}_p{}".format(a, f, p)) # Old line, works
              dataOfFuction.append(tempData)
              break
    # Each list from solutions contains the data for each function
    solutions.append(dataOfFuction)
  return solutions

# Get style (css) for tables
def getTableStyle():
  th_props = [
      ('font-size', '18px'),
      ('font-family', 'monospace'),
      ('text-align', 'center'),
      # ('background-color', '#f7f7f9'),
      ('text-align', 'center')
  ]

  td_props = [
    ('font-size', '16px'),
    ('font-family', 'monospace'),
    ('padding', '10px'),
    ('padding-top', '25px'),
    ('text-align', 'center')
  ]

  title_props = [
    ('font-size', '24px'),
    ('font-family', 'monospace'),
    ('text-align', 'center'),
    ('font-weight', 'bold'),
    ("margin-bottom", "30px"),
    ("margin-top", "30px")
  ]

  styles = [
    dict(selector="th", props=th_props),
    dict(selector="td", props=td_props),
    dict(selector="caption", props=title_props)
  ]
  return styles

def getKrampserResults(grouped, index, case):
  # Maxfe = 15k continuous
  if index == 0: # t10
    if case == "continuous":
      tableTitle = "10-bar truss problem, NEvals: 15000 | DUVDE e APM: NEvals: 280000"
      grouped.loc[4] = ["SMDE k=2*", 5060.87, 5060.92, 5061.98, 3.93e+00, 5076.70] # Append row at dataFrame
      grouped.loc[5] = ["DUVDE*", 5060.85, None, 5067.18, 7.94e+00, 5076.66] # Append row at dataFrame
      grouped.loc[6] = ["APM*", 5069.08, None, 5091.43, None, 5117.39] # Append row at dataFrame
    elif case == "discrete":
      tableTitle = "10-bar truss problem, NEvals: 15000 | DUVDE NEvals: 24000 e APM NEvals: 90000"
      grouped.loc[4] = ["SMDE k=2*", 5490.74, 5490.74, 5495.99, 1.13e+01, 5529.30] # Append row at dataFrame
      grouped.loc[5] = ["DUVDE*", 5562.35, None, 5564.90, 0.6, 5565.04] # Append row at dataFrame  24k
      grouped.loc[6] = ["APM*", 5490.74, None, 5545.48, None, 5567.84] # Append row at dataFrame 90k
    else:
      sys.exit("Case not defined.")
  elif index == 1: # t25
    if case == "continuous":
      tableTitle = "25-bar truss problem, NEvals: 15000"
      grouped.loc[4] = ["SMDE k=2*", 484.06, 484.07, 484.07, 0.0107, 484.10] # Append row at dataFrame
    elif case == "discrete":
      tableTitle = "25-bar truss problem, NEvals: 15000 | DUVDE e APM: NEvals: 20000"
      grouped.loc[4] = ["SMDE k=2*", 484.85, 485.05, 485.44, 0.693, 487.13] # Append row at dataFrame
      grouped.loc[5] = ["DUVDE*", 485.90, None, 498.44, 7.66e+00, 507.77] # Append row at dataFrame 20k
      grouped.loc[6] = ["APM*", 485.85, None, 485.97, None, 490.74] # Append row at dataFrame 20k
  elif index == 2: # t60
    if case == "continuous":
      tableTitle = "60-bar truss problem, NEvals: 15000 | DUVDE NEvals:150000 e APM NEvals: 800000"
      grouped.loc[4] = ["SMDE k=2*", 308.94, 309.42, 309.49, 0.464, 311.21] # Append row at dataFrame
      grouped.loc[5] = ["DUVDE", 309.44, None, 311.54, 1.46e+00, 314.70] # Append row at dataFrame
      grouped.loc[6] = ["APM", 311.87, None, 333.01, None, 384.19] # Append row at dataFrame
    elif case == "discrete":
      tableTitle = "60-bar truss problem, NEvals: 15000"
      grouped.loc[4] = ["SMDE k=2*", 312.73, 314.20, 315.12, 3.98e+00, 335.88] # Append row at dataFrame
    else:
      sys.exit("Case not defined.")
  elif index == 3: # t72
    if case == "continuous":
      tableTitle = "72-bar truss problem, NEvals: 15000 | DUVDE e APM: NEvals: 35000"
      grouped.loc[4] = ["SMDE k=2*", 379.62, 379.63, 379.65, 0.0341, 379.73] # Append row at dataFrame
      grouped.loc[5] = ["DUVDE", 379.66, None, 380.42, 0.572, 381.37] # Append row at dataFrame
      grouped.loc[6] = ["APM", 387.04, None, 402.59, None, 432.95] # Append row at dataFrame
    elif case == "discrete":
      tableTitle = "72-bar truss problem, NEvals: 15000"
      grouped.loc[4] = ["SMDE k=2*", 385.54, 386.81, 386.91, 1.05e+00, 389.21] # Append row at dataFrame
    else:
      sys.exit("Case not defined.")
  elif index == 4: # t942
    if case == "continuous":
      tableTitle = "942-bar truss problem, NEvals: 15000"
      grouped.loc[4] = ["SMDE k=2*", 149932.00, 171218.50, 174369.63, 1.72e+04, 230139.00] # Append row at dataFrame
    elif case == "discrete":
      tableTitle = "942-bar truss problem, NEvals: 15000"
      grouped.loc[4] = ["SMDE k=2*", 153010.00, 181899.00, 181127.23, 1.73e+04, 216514.00] # Append row at dataFrame
    else:
      sys.exit("Case not defined.")
  else:
    sys.exit("Index not defined.")
  return grouped, tableTitle

def getExtraResults(np2dArr, rowsTitles, index, case):
  # Maxfe = 15k continuous
  if index == 0: # t10
    if case == "continuous":
      # tableTitle = "10-bar truss problem, NEvals: 15000 | DUVDE e APM: NEvals: 280000"
      # Define algorithms and results
      smde = ["SMDE k=2*", 5060.87, 5060.92, 5061.98, 3.93e+00, 5076.70] 
      duvde = ["DUVDE*", 5060.85, float("NaN"), 5067.18, 7.94e+00, 5076.66] 
      apm = ["APM*", 5069.08, float("NaN"), 5091.43, float("NaN"), 5117.39]

      # Append name of algorihtms
      rowsTitles.append(smde[0])
      rowsTitles.append(duvde[0])
      rowsTitles.append(apm[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [duvde[1:]], axis=0)
      np2dArr = np.append(np2dArr, [apm[1:]], axis=0)
    elif case == "discrete":
      smde = ["SMDE k=2*", 5490.74, 5490.74, 5495.99, 1.13e+01, 5529.30]
      duvde = ["DUVDE*", 5562.35, float("NaN"), 5564.90, 0.6, 5565.04]
      apm = ["APM*", 5490.74, float("NaN"), 5545.48, float("NaN"), 5567.84]

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
  elif index == 1: # t25
    if case == "continuous":
      # tableTitle = "25-bar truss problem, NEvals: 15000"
      smde= ["SMDE k=2*", 484.06, 484.07, 484.07, 0.0107, 484.10]
      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0) # Appends values of given list (remove first index, which is the string)
    elif case == "discrete":
      # tableTitle = "25-bar truss problem, NEvals: 15000 | DUVDE e APM: NEvals: 20000"
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
  elif index == 2: # t60
    if case == "continuous":
      # tableTitle = "60-bar truss problem, NEvals: 15000 | DUVDE NEvals:150000 e APM NEvals: 800000"
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
    elif case == "discrete":
      # tableTitle = "60-bar truss problem, NEvals: 15000"
      smde = ["SMDE k=2*", 312.73, 314.20, 315.12, 3.98e+00, 335.88]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  elif index == 3: # t72
    if case == "continuous":
      # tableTitle = "72-bar truss problem, NEvals: 15000 | DUVDE e APM: NEvals: 35000"
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
    elif case == "discrete":
      # tableTitle = "72-bar truss problem, NEvals: 15000"
      smde = ["SMDE k=2*", 385.54, 386.81, 386.91, 1.05e+00, 389.21]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  elif index == 4: # t942
    if case == "continuous":
      # tableTitle = "942-bar truss problem, NEvals: 15000"
      smde = ["SMDE k=2*", 149932.00, 171218.50, 174369.63, 1.72e+04, 230139.00]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    elif case == "discrete":
      tableTitle = "942-bar truss problem, NEvals: 15000"
      smde = ["SMDE k=2*", 153010.00, 181899.00, 181127.23, 1.73e+04, 216514.00]

      # Append name of algorihtms
      rowsTitles.append(smde[0])

      # Append reults on np2darr
      np2dArr = np.append(np2dArr, [smde[1:]], axis=0)
    else:
      sys.exit("Case not defined.")
  else:
    sys.exit("Index not defined.")
  return np2dArr, rowsTitles

# Generates in html
def makeAnalysisHtml(solutions):
  # Read solutions and store on dataframe
  case="discrete"
  for index, solution in enumerate(solutions):
    df = pd.DataFrame.from_records(solution)
    titles = ["Objective Function", "ViolationSum/Fitness"]
    # Define name for columns of the project variables
    for i in range(len(df.columns)-4):
      titles.append("ProjVar{}".format(i))
    # Define name for the last 2 columns
    titles.append("CPU Time(s)")
    titles.append("Algorithm")
    # Name columns
    df.columns = titles

    # Drop columns that contains word | Could be done with df.filter(regex)
    df.drop([col for col in df.columns if "ProjVar" in col], axis=1, inplace=True)
    df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
    df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)

    # pd.set_option("display.float_format","{:e}".format) # Use scientific notation as format
    pd.set_option("display.float_format","{:.2f}".format)
    # Group by algorithm 
    grouped = df.groupby(["Algorithm"])
    grouped = grouped.agg(["min", "median", "mean", np.std, "max"]).reset_index(level="Algorithm") # Need to reset index for append new column
    # grouped = grouped.agg(["min", "median", "mean", "std", "max"]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works

    # print(grouped)
    # grouped.loc[4] = ["Alg1", 5070.00, 15, 5060, 5080] # Append row at dataFrame
    # Maxfe = 15k continuous
    # if(index == 0): # t10
    #   grouped.loc[4] = ["SMDE k=2*", 5060.87, 5060.92, 5061.98, 3.93e+00, 5076.70] # Append row at dataFrame
    #   grouped.loc[5] = ["DUVDE*", 5060.85, None, 5067.18, 7.94e+00, 5076.66] # Append row at dataFrame
    #   grouped.loc[6] = ["APM*", 5069.08, None, 5091.43, None, 5117.39] # Append row at dataFrame
    # elif(index == 1): # t25
    #   grouped.loc[4] = ["SMDE k=2*", 484.06, 484.07, 484.07, 0.0107, 484.10] # Append row at dataFrame
    # elif(index == 2): # t60
    #   grouped.loc[4] = ["SMDE k=2*", 308.94, 309.42, 309.49, 0.464, 311.21] # Append row at dataFrame
    #   grouped.loc[5] = ["DUVDE", 309.44, None, 311.54, 1.46e+00, 314.70] # Append row at dataFrame
    #   grouped.loc[6] = ["APM", 311.87, None, 333.01, None, 384.19] # Append row at dataFrame
    # elif(index == 3): # t72
    #   grouped.loc[4] = ["SMDE k=2*", 379.62, 379.63, 379.65, 0.0341, 379.73] # Append row at dataFrame
    #   grouped.loc[5] = ["DUVDE", 379.66, None, 380.42, 0.572, 381.37] # Append row at dataFrame
    #   grouped.loc[6] = ["APM", 387.04, None, 402.59, None, 432.95] # Append row at dataFrame
    # elif(index == 4): # t942
    #   grouped.loc[4] = ["SMDE k=2*", 149932.00, 171218.50, 174369.63, 1.72e+04, 230139.00] # Append row at dataFrame
    # else:
    #   sys.exit("Index not defined.")


    # # Maxfe = 15k discrete
    # if(index == 0): # t10
    #   grouped.loc[4] = ["SMDE k=2*", 5490.74, 5490.74, 5495.99, 1.13e+01, 5529.30] # Append row at dataFrame
    #   grouped.loc[5] = ["DUVDE*", 5562.35, None, 5564.90, 0.6, 5565.04] # Append row at dataFrame  24k
    #   grouped.loc[6] = ["APM*", 5490.74, None, 5545.48, None, 5567.84] # Append row at dataFrame 90k
    # elif(index == 1): # t25
    #   grouped.loc[4] = ["SMDE k=2*", 484.85, 485.05, 485.44, 0.693, 487.13] # Append row at dataFrame
    #   grouped.loc[5] = ["DUVDE*", 485.90, None, 498.44, 7.66e+00, 507.77] # Append row at dataFrame 20k
    #   grouped.loc[6] = ["APM*", 485.85, None, 485.97, None, 490.74] # Append row at dataFrame 20k
    # elif(index == 2): # t60
    #   grouped.loc[4] = ["SMDE k=2*", 312.73, 314.20, 315.12, 3.98e+00, 335.88] # Append row at dataFrame
    # elif(index == 3): # t72
    #   grouped.loc[4] = ["SMDE k=2*", 385.54, 386.81, 386.91, 1.05e+00, 389.21] # Append row at dataFrame
    # elif(index == 4): # t942
    #   grouped.loc[4] = ["SMDE k=2*", 153010.00, 181899.00, 181127.23, 1.73e+04, 216514.00] # Append row at dataFrame
    # else:
    #   sys.exit("Index not defined.")

    grouped, tableTitle = getKrampserResults(grouped, index, case)

    # grouped = grouped.drop("Algorithm", axis=1) # Drop algorithm columns

    # grouped.index.names = ["Algorithm"] # Adds name on index column
    # grouped.index=(["a", "b", "c", "d", "e"]) # Rename indexes from

    # print(grouped)

    # Gets styles for formatting table and apply on df
    
    styles = getTableStyle()
    t = (grouped.style
          .set_table_styles(styles)
          .apply(highlight_min)
          # .highlight_min()
          .set_caption(tableTitle)
          # .set_precision(2)
          )

    # Generates html+css table from dataframe
    
    print(t.render()) 
    # break
    
    # df.to_csv('tableCsv.csv') # Exports dataframe to csv file
    # df.to_latex('tableLatex.txt') #  Export dataframe to simple latex table

# Generates statistical measures and generate .tex table
def makeAnalysis(solutions):
  # Read solutions and store on dataframe
  case="discrete"
  for key, solution in enumerate(solutions):
    df = pd.DataFrame.from_records(solution)
    titles = ["Objective Function", "ViolationSum/Fitness"]
    # Define name for columns of the project variables
    for i in range(len(df.columns)-4):
      titles.append("ProjVar{}".format(i))

    # Define name for the last 2 columns
    titles.append("CPU Time(s)")
    titles.append("Algorithm")

    # Name columns
    df.columns = titles

    # Drop columns that contains word | Could be done with df.filter(regex)
    df.drop([col for col in df.columns if "ProjVar" in col], axis=1, inplace=True)
    df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
    df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)

    # Group by algorithm and calcualtes statistical measures
    grouped = df.groupby(["Algorithm"])
    grouped = grouped.agg(["min", "median", "mean", "std", "max"]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works

    # Generates 2d np array from pandas grouped df
    np2dArr = grouped.to_numpy()
    # Get row and columns titles
    rowsTitles = list(grouped.index.values) 
    columnsTitles = list(grouped.columns.levels[1])
    # Get extra results (appending other algorithms on 2d np arr)
    np2dArr, rowsTitles = getExtraResults(np2dArr, rowsTitles, key, case)
    # Generates .text table
    np2dArrToLatex(np2dArr, columnsTitles, rowsTitles)

    # break

def np2dArrToLatex(np2dArr, columnsTitles, rowsTitles):
  # List containing minimum value of each column
  minOfEachColumn = np2dArr.min(axis=0) 

  # Define print format, table aligmnent and round rpecision
  precFormat = ".2f" # Precision format for priting table. Cand be '.xf' or 'e' for scientific notation
  texTableAlign = "r " # Align on right (could be l, c, r and so on)
  tabAlign = "{{@{{}} l | {} @{{}}}}".format(texTableAlign*len(columnsTitles)) 
  precRound = 2 # Precision round parameter for highlightning minimum value from column

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
  print("\\textbf{{{}}}: {}".format("Caso", "Continuo"))
  print("\\end{table}")
  # Begin of table structure .tex
  print()

if __name__ == '__main__':
  solutions = readResults("last")
  makeAnalysis(solutions)



# *---*---*---*---*---*---*---*---*---* Pandas Personal Helper ---*---*---*---*---*---*---*---*---*---*---*
# print(grouped.index.names) # Return indexes names, in this case, only "Algorithm" is returned
# print(grouped.index) # Each index is each algorithm "CMAES+APM", "CMAES+DEB" ...
# print(grouped.index.values) # Returns name of all algorithms

# print(grouped.columns.levels) # Return a list of list containing all column levels, from top to bottom.
# print(list(grouped.columns.levels[1].values)) # Return all names of columns given a index
import numpy as np
import pandas as pd
import os
from pathlib import Path

def readResults():
  algorithms = ["DE", "CMAES"]
  # functions = [10, 25, 60, 72, 942]
  # functions = [10]
  functions = [10, 25, 60, 72]
  constraintHandlingMethods = [1, 2]
  # seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  # , 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
  seeds = [1, 2, 3, 4, 5]

  basePath =  Path(__file__).parent
  # filePath = (basePath / "../results/t10").resolve()

  solutions = []
  # solutions = [[dadosf10], [dadosf25]]
  # for f in (functions): # Functions
  for f in (functions): # Functions
    tempFunc = []
    for a in (algorithms): # Algorithms
      for p in (constraintHandlingMethods): # Constraint handling methods
        for s in (seeds): # Seeds
          tempData = []
          filePath = (basePath / "../results/t{}/{}_f{}_p{}_s{}.dat".format(f, a, f, p, s)).resolve()
          file = open(filePath)
          while True:
            buffer = file.readline()
            if "Hall of Fame" in buffer: # Find hall of fame line
              buffer = file.readline() # Read one more
              buffer = buffer.split(" ")
              for item in buffer:
                tempData.append(float(item))
              
              # print(type(solutions))
              # print(solutions)
            elif "CPU time used" in buffer:
              buffer = file.readline()
            #   # Inserts time at the end of each solution (after project variables)
            #   solutions[len(solutions)-1].append(cpuTimeUsed)
              # solutions[len(solutions)-1]+=buffer
              tempData.append(float(buffer))
              tempData.append("{}_f{}_p{}".format(a, f, p))
              # solutions.append(tempData)
              tempFunc.append(tempData)
              break
    solutions.append(tempFunc)
  return solutions

# Puts bold on the maximum value of column
def highlight_max(s):
  is_max = s == s.max()
  return ['font-weight: bold' if v else '' for v in is_max]

# Puts bold on the minimum value of column
def highlight_min(s):
  is_max = s == s.min()
  return ['font-weight: bold' if v else '' for v in is_max]

# Get style (css) for tables
def getTableStyle():
  th_props = [
      ('font-size', '18px'),
      ('font-family', 'monospace'),
      ('text-align', 'center'),
      ('background-color', '#f7f7f9'),
  ]

  td_props = [
    ('font-size', '16px'),
    ('font-family', 'monospace'),
    ('padding', '10px'),
    ('padding-top', '25px')
  ]

  styles = [
    dict(selector="th", props=th_props),
    dict(selector="td", props=td_props)
  ]
  return styles

def makeAnalysisNew(solutions):
  # Read solutions and store on dataframe
  for solution in solutions:
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

    # Rename all rows from column "Algorithm" that contains searched string (more readibility)
    # df["Algorithm"] = df["Algorithm"].str.replace("DE_f10_p1", "DE + Deb")
    # df["Algorithm"] = df["Algorithm"].str.replace("DE_f10_p2", "DE + APM")
    # df["Algorithm"] = df["Algorithm"].str.replace("CMAES_f10_p1", "CMAES + Deb")
    # df["Algorithm"] = df["Algorithm"].str.replace("CMAES_f10_p2", "CMAES + APM")
    # Drop columns that contains word | Could be done with df.filter(regex)
    df.drop([col for col in df.columns if "ProjVar" in col], axis=1, inplace=True)
    df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
    df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)

    # pd.set_option("display.float_format","{:e}".format) # Use scientific notation as format
    # Group by algorithm 
    df = df.groupby(["Algorithm"])
    
    # Aggregate and realize the math
    df = df.agg(["mean", "std", "min", "max"]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works



    # Gets styles for formatting table and apply on df
    styles = getTableStyle()
    t = (df.style
          .set_table_styles(styles)
          .apply(highlight_min))



    # Generates html+css table from dataframe
    print(t.render()) 
    
    # df.to_csv('tableCsv.csv') # Exports dataframe to csv file
    # df.to_latex('tableLatex.txt') #  Export dataframe to simple latex table

def makeAnalysis(solutions):
  # Read solutions and store on dataframe
  df = pd.DataFrame.from_records(solutions)
  titles = ["Objective Function", "ViolationSum/Fitness"]
  # Define name for columns of the project variables
  for i in range(len(df.columns)-4):
    titles.append("ProjVar{}".format(i))
  # Define name for the last 2 columns
  titles.append("CPU Time(s)")
  titles.append("Algorithm")
  # Name columns
  df.columns = titles

  # Rename all rows from column "Algorithm" that contains searched string (more readibility)
  df["Algorithm"] = df["Algorithm"].str.replace("DE_f10_p1", "DE + Deb")
  df["Algorithm"] = df["Algorithm"].str.replace("DE_f10_p2", "DE + APM")
  df["Algorithm"] = df["Algorithm"].str.replace("CMAES_f10_p1", "CMAES + Deb")
  df["Algorithm"] = df["Algorithm"].str.replace("CMAES_f10_p2", "CMAES + APM")
  # Drop columns that contains word | Could be done with df.filter(regex)
  df.drop([col for col in df.columns if "ProjVar" in col], axis=1, inplace=True)
  df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
  df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)

  # pd.set_option("display.float_format","{:e}".format) # Use scientific notation as format
  # Group by algorithm 
  df = df.groupby(["Algorithm"])
  
  # Aggregate and realize the math
  df = df.agg(["mean", "std", "min", "max"]) # df = df.agg([np.mean, np.std, np.min, np.max]) also works



  # Gets styles for formatting table and apply on df
  styles = getTableStyle()
  t = (df.style
        .set_table_styles(styles)
        .apply(highlight_min))



  # Generates html+css table from dataframe
  print(t.render()) 
  
  # df.to_csv('tableCsv.csv') # Exports dataframe to csv file
  # df.to_latex('tableLatex.txt') #  Export dataframe to simple latex table


if __name__ == '__main__':
  solutions = readResults()
  # print(solutions)
  # print(len(solutions))
  # print("solutions 0")
  # print(solutions[0])
  # print(solutions[1])
  makeAnalysisNew(solutions)
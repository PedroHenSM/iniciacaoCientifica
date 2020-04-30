import numpy as np
import pandas as pd
import os
from pathlib import Path

def readResults():
  algorithms = ["DE", "CMAES"]
  # functions = [10, 25, 60, 72, 942]
  # functions = [10, 25, 60, 72]
  functions = [10]
  constraintHandlingMethods = [1, 2]
  # seeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
  # , 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
  seeds = [1, 2, 3, 4, 5]

  basePath =  Path(__file__).parent
  # filePath = (basePath / "../results/t10").resolve()

  solutions = []
  for f in (functions): # Functions
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
              solutions.append(tempData)
              break
  return solutions

def makeAnalysis(solutions):
  df = pd.DataFrame.from_records(solutions)
  titles = ["Objective Function", "ViolationSum/Fitness"]
  for i in range(len(df.columns)-4):
    titles.append("ProjVar{}".format(i))
  titles.append("CPU Time(s)")
  titles.append("Alg")
  df.columns = titles
  print(df)
  # print(df.describe())
  # print(df.iloc[4]) # Selects 4 row of the df
  toDrop = "ProjVar"
  # Drop columns that contains word "toDrop" | Could be done with df.filter(regex)
  df.drop([col for col in df.columns if toDrop in col], axis=1, inplace=True)
  df.drop([col for col in df.columns if "ViolationSum" in col], axis=1,inplace=True)
  df.drop([col for col in df.columns if "CPU Time" in col], axis=1,inplace=True)
  print("After drop columns")
  print(df)
  df = df.groupby(["Alg"])
  print("After grouping")
  print(df)
  # df = df.describe()
  # df = df.mean()
  # print("Describe")
  # print(df)
  # df = df.describe()
  # df.describe().unstack()[['count','max']]
  # df = df.describe().loc[['count','max']]
  df1 = df.agg([np.mean, np.std, np.min, np.max]) # Works
  df2 = df.agg(["mean", "std", "min", "max"]) # Works

  print("Agg")
  print(df1)
  print(df2)
  # print(a)

  # print(a)

if __name__ == '__main__':
  solutions = readResults()
  # print(solutions)
  # print("solutions 0")
  # print(solutions[0])
  # print(solutions[1])
  makeAnalysis(solutions)
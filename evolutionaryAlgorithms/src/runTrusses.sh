#!/bin/bash

totalAlgorithms=2 # 2
algorithms=(DE CMAES)

totalFunctions=4 # 5
functions=(110 125 160 172 1942)

totalSeeds=5 # 30
seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)

totalConstraintHandlingMethods=2 # 2 | 1 Deb and 2 APM
constraintHandlings=(DEB APM)

totalPopulationSize=1 # ParentsSize
populations=(50 100)

totalOffspringsSize=1 # offspringsSize
offsprings=(50 100 200 400)

totalFe=1 # Max functions evaluations
functionEvaluations=(20000)

f=0
while(($f<$totalFunctions))
do
  l=0
  while(($l<$totalOffspringsSize))
  do
    u=0
    while(($u<$totalPopulationSize))
    do
      fe=0
      while(($fe<$totalFe))
      do
        p=0
        while(($p<$totalConstraintHandlingMethods))
        do
          s=0
          while(($s<$totalSeeds))
          do
            a=0
            while(($a<$totalAlgorithms))
            do
              echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} "
              python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} > \
              ../results/t${functions[f]:1}/${algorithms[a]}_f${functions[f]:1}_p${constraintHandlings[p]}_s${seeds[s]}.dat
              a=$((a+1))
            done
            s=$((s+1))
          done
          p=$((p+1))
        done
        fe=$((fe+1))
      done
      u=$((u+1))
    done
    l=$((l+1))
  done
  f=$((f+1))
done







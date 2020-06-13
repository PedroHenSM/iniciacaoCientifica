#!/bin/bash

algorithms=(DE CMAES) # Algorithms
totalAlgorithms=${#algorithms[@]} 

functions=(21 22 23 24 25 110 125 160 172 1942) # Problems
totalFunctions=${#functions[@]} # 5

seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30) # Seeds
totalSeeds=${#seeds[@]}

constraintHandlings=(DEB APM) # Constraint handling methods
totalConstraintHandlingMethods=${#constraintHandlings[@]}

populations=(50) # Populations (parents) size
totalPopulationSize=${#populations[@]}

offsprings=(50) # Populations (offsprings) size
totalOffspringsSize=${#offsprings[@]}

weights=(Linear Superlinear Equal) # CMA-ES weights parameters
totalWeights=${#weights[@]}

functionEvaluations=(15000) # Trusses
totalFe=${#functionEvaluations[@]} # Max functions evaluations

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
          w=0
          while(($w<$totalWeights))
          do
            s=0
            while(($s<$totalSeeds))
            do
              a=0
              while(($a<$totalAlgorithms))
              do
                # Executing CMAES, should run with all three weights parameters
                if [ ${algorithms[a]} = CMAES ]
                then
                  if [ ${functions[f]} -eq 21 -o ${functions[f]} -eq 22 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 36000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 36000 -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                  # Welded beam has a max evaluation of 320k
                  elif [ ${functions[f]} -eq 23 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 320000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 320000 -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                  # Welded beam has a max evaluation of 80k
                  elif [ ${functions[f]} -eq 24 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 80000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 80000 -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                  # Cantilever beam has a max evaluation of 35k
                  elif [ ${functions[f]} -eq 25 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 35000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 35000 -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                  # Other functions read max evaluation from functionEvaluations array
                  else 
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                  fi
                # Only executes DE when using the first weights(this parameters doesnt change DE results, no need to run with 3 different ones)
                elif [ ${algorithms[a]} = DE -a $w -eq 0 ]
                then
                  if [ ${functions[f]} -eq 21 -o ${functions[f]} -eq 22 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 36000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 36000 > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                  # Welded beam has a max evaluation of 320k
                  elif [ ${functions[f]} -eq 23 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 320000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 320000 > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                  # Welded beam has a max evaluation of 80k
                  elif [ ${functions[f]} -eq 24 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 80000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 80000 > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                  # Cantilever beam has a max evaluation of 35k
                  elif [ ${functions[f]} -eq 25 ];
                  then
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 35000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 35000 > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                  # Other functions read max evaluation from functionEvaluations array
                  else 
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                  fi
                fi
                a=$((a+1))
              done
              s=$((s+1))
            done
            w=$((w+1))
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







#!/bin/bash

# Function 21 => Tension compression spring
# Function 22 => Speed reducer
# Function 23 => Welded beam
# Function 24 => Pressure vessel
# Function 25 => Cantilever beam

start=$(date +"%c") # Get curent time

algorithms=(DE CMAES) # Algorithms
totalAlgorithms=${#algorithms[@]} 


functions=(21 22 23 24 25) # Problems
totalFunctions=${#functions[@]}

seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30) # Seeds
totalSeeds=${#seeds[@]}

problemsCase=(continuous discrete)
totalProblemsCase=${#problemsCase[@]}

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
                  # Tension/compression spring and Speed reducer
                  if [ ${functions[f]} -eq 21 -o ${functions[f]} -eq 22 ]
                  then
                    maxEval=36000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                  # Welded beam has a max evaluation of 320k
                  elif [ ${functions[f]} -eq 23 ]
                  then
                    maxEval=320000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                  # Pressure vesel has a max evaluation of 80k
                  elif [ ${functions[f]} -eq 24 ]
                  then
                    maxEval=80000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                  # Cantilever beam has a max evaluation of 35k
                  elif [ ${functions[f]} -eq 25 ]
                  then
                    maxEval=35000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                  fi
                # Only executes DE when using the first weights(these parameters doesnt change DE results, no need to run with 3 different ones)
                elif [ ${algorithms[a]} = DE -a $w -eq 0 ]
                then
                  # Tension/compression spring and Speed reducer
                  if [ ${functions[f]} -eq 21 -o ${functions[f]} -eq 22 ]
                  then
                    maxEval=36000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                  # Welded beam has a max evaluation of 320k
                  elif [ ${functions[f]} -eq 23 ]
                  then
                    maxEval=320000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                  # Pressure vesel has a max evaluation of 80k
                  elif [ ${functions[f]} -eq 24 ]
                  then
                    maxEval=80000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                  # Cantilever beam has a max evaluation of 35k
                  elif [ ${functions[f]} -eq 25 ]
                  then
                    maxEval=35000
                    echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                    python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                    ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
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

end=$(date +"%c") # Get curent time
echo "Started execution at: $start"
echo "Finished execution at: $end "







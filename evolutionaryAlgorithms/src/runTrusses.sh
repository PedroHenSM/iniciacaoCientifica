#!/bin/bash

# 10 bar truss - Continuous: 280k | Discrete: 240k and 90k
# 25 bar truss - Continuous: 240k | Discrete: 20k
# 60 bar truss - Continuous: 12k (Rafael), 150k e 800k (Krempser)| Discrete:
# 72 bar truss - Continuous: 35k | Discrete:
# 942 bar truss -  Continuous: | Discrete:

start=$(date +"%c") # Get curent time

algorithms=(DE CMAES) # Algorithms
totalAlgorithms=${#algorithms[@]} 

# functions=(110 125 160 172 1942) # Problems
functions=(110 ) # Problems
totalFunctions=${#functions[@]}

seeds=(1 2) # Seeds
# seeds=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30) # Seeds
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
  c=0
  while(($c<$totalProblemsCase))
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
                    # 10 bar truss
                    if [ ${functions[f]} -eq 110 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=280000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        maxEval=90000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 25 bar truss
                    elif [ ${functions[f]} -eq 125 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=240000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        maxEval=20000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 60 bar truss
                    elif [ ${functions[f]} -eq 160 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=12000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 72 bar truss
                    elif [ ${functions[f]} -eq 172 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=35000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 942 bar truss
                    elif [ ${functions[f]} -eq 1942 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} -w ${weights[w]} -c ${problemsCase[c]} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    fi
                  # Only executes DE when using the first weights(these parameters doesnt change DE results, no need to run with 3 different ones)
                  elif [ ${algorithms[a]} = DE -a $w -eq 0 ]
                  then
                    # 10 bar truss
                    if [ ${functions[f]} -eq 110 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=280000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        maxEval=90000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 25 bar truss
                    elif [ ${functions[f]} -eq 125 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=240000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        maxEval=20000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 60 bar truss
                    elif [ ${functions[f]} -eq 160 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=12000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 72 bar truss
                    elif [ ${functions[f]} -eq 172 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        maxEval=35000
                        echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      fi
                    # 942 bar truss
                    elif [ ${functions[f]} -eq 1942 ]
                    then
                      # Continuous case
                      if [ ${problemsCase[c]::1} = 'c' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      # Discrete case
                      elif [ ${problemsCase[c]::1} = 'd' ]
                      then
                        : # Does nothing
                        # maxEval=280000
                        # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${maxEval} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                        # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${maxEval} > \
                        # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_m${maxEval}_s${seeds[s]}.dat
                      fi
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
    c=$((c+1))
  done
  f=$((f+1))
done

end=$(date +"%c") # Get curent time
echo "Started execution at: $start"
echo "Finished execution at: $end "







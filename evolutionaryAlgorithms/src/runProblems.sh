#!/bin/bash

# 10 bar truss - Continuous: 280k | Discrete: 90k and 24k
# 25 bar truss - Continuous: 240k | Discrete: 20k
# 60 bar truss - Continuous: 12k (Rafael), 150k e 800k (Krempser)| Discrete:
# 72 bar truss - Continuous: 35k | Discrete:
# 942 bar truss - Continuous: | Discrete:
# Function 21 Tension compression spring
# Function 22 Speed reducer
# Function 23 Welded beam
# Function 24 Pressure vessel
# Function 25 Cantilever beam

start=$(date +"%T") # Get curent time

algorithms=(DE CMAES) # Algorithms
totalAlgorithms=${#algorithms[@]} 

functions=(21 22 23 24 25 110 125 160 172 1942) # Problems
# functions=(21 22 23 24 25 110 125 160 172) # Problems
# functions=(1942) # Problems
# functions=(21) # Problems
# functions=(21 22 23 24 25) # Problems
# functions=(110 125 160 172) # Problems
totalFunctions=${#functions[@]}

# seeds=(1 2) # Seeds
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

# weights=(Superlinear) # CMA-ES weights parameters
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
                    # Only executes classical engineering problems when using the first problem type (this parameters doesnt change these problems results, no need to run with 2 different ones)
                    if [ ${functions[f]} -eq 21 -a $c -eq 0 -o ${functions[f]} -eq 22 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 36000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 36000 -w ${weights[w]} > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    # Welded beam has a max evaluation of 320k
                    elif [ ${functions[f]} -eq 23 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 320000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 320000 -w ${weights[w]} > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    # Welded beam has a max evaluation of 80k
                    elif [ ${functions[f]} -eq 24 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 80000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 80000 -w ${weights[w]} > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    # Cantilever beam has a max evaluation of 35k
                    elif [ ${functions[f]} -eq 25 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 35000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 35000 -w ${weights[w]} > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    # Other functions read max evaluation from functionEvaluations array
                    elif [ ${functions[f]::1} -eq 1 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -w ${weights[w]} -c ${problemsCase[c]} > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    else 
                      echo "Run nothing"
                      # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                      # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -w ${weights[w]} -c ${problemsCase[c]} > \
                      # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    fi
                  # Only executes DE when using the first weights(this parameters doesnt change DE results, no need to run with 3 different ones)
                  elif [ ${algorithms[a]} = DE -a $w -eq 0 ]
                  then
                    if [ ${functions[f]} -eq 21 -a $c -eq 0 -o ${functions[f]} -eq 22 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 36000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 36000 > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                    # Welded beam has a max evaluation of 320k
                    elif [ ${functions[f]} -eq 23 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 320000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 320000 > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                    # Welded beam has a max evaluation of 80k
                    elif [ ${functions[f]} -eq 24 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 80000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 80000 > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                    # Cantilever beam has a max evaluation of 35k
                    elif [ ${functions[f]} -eq 25 -a $c -eq 0 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: 35000 parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m 35000 > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_p${constraintHandlings[p]}_s${seeds[s]}.dat
                    # Other functions read max evaluation from functionEvaluations array
                    elif [ ${functions[f]::1} -eq 1 ];
                    then
                      echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                      python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -w ${weights[w]} -c ${problemsCase[c]} > \
                      ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_w${weights[w]}_s${seeds[s]}.dat
                    else 
                      echo "Run nothing"
                      # echo "Executing algorithm: ${algorithms[a]} seed: ${seeds[s]} constraint handling method: ${constraintHandlings[p]} maxFe: ${functionEvaluations[fe]} parentsSize: ${populations[u]} offspringsSize: ${offsprings[l]} function: ${functions[f]} weights: ${weights[w]} case: ${problemsCase[c]}"
                      # python3 main.py -a ${algorithms[a]} -f ${functions[f]} -s ${seeds[s]} -p ${constraintHandlings[p]} -u ${populations[u]} -l ${offsprings[l]} -m ${functionEvaluations[fe]} -c ${problemsCase[c]} > \
                      # ../results/functions/f${functions[f]}/${algorithms[a]}_f${functions[f]}_c${problemsCase[c]::1}_p${constraintHandlings[p]}_s${seeds[s]}.dat
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

end=$(date +"%T") # Get curent time
echo "Started execution at: $start"
echo "Finished execution at: $end "







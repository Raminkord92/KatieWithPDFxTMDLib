#! /bin/bash

set -e

#Begin process list
Nproc=14
array[0]="echo A 1"
array[1]="echo A 2"
array[2]="echo A 3"
array[3]="echo B 1"
array[4]="echo B 2"
array[5]="echo B 3"
array[6]="echo C 1"
array[7]="echo C 2"
array[8]="echo C 3"
array[9]="echo D 1"
array[10]="echo D 2"
array[11]="echo D 3"
array[12]="echo E 1"
array[13]="echo E 2"
#End process list

for ((i=0; i<$Nproc; i++)); do
  proc[$i]=${array[$i]}
done

Ncpu=$Nproc

for var in "$@"
do
    lhs=`echo $var | sed s~=.*$~~`
    rhs=`echo $var | sed s~^$lhs=~~`
    case $lhs in
    "help"|"-h"|"-help"|"--help"|"-i"|"-info"|"--info")
        echo "==============================================="
        echo " To run all optimization processes at once:"
        echo " $ ./optimize.sh"
        echo " To run 4 optimization processes at a time:"
        echo " $ ./optimize.sh Ncpu=4"
        echo " Exactly the same is achieved with Nparallel=4"
        echo " To run optimization process 3:"
        echo " $ ./optimize.sh proc=3"
        echo " To run optimization process 3 and 12:"
        echo " $ ./optimize.sh proc=3,12"
        echo " You should now understand the following:"
        echo " $ ./optimize.sh proc=3,12,4,11,2 Ncpu=4"
        echo " You can monitor the progress with"
        echo " $ tail -f proc*/output"
        echo " You can kill all processes with"
        echo " $ pkill -f main"
        echo "==============================================="
        exit 0
        ;;
    "Ncpu")
        Ncpu=$rhs
        ;;
    "Nparallel")
        Ncpu=$rhs
        ;;
    "proc")
        list=(`echo $rhs | sed -e's/,/ /g'`)
        Nproc=${#list[@]}
        for ((i=0; i<$Nproc; i++)); do
          proc[$i]=${array[${list[$i]}-1]}
        done
        ;;
    esac
done

if [ "$#" -eq 0 ]; then
  Nsugg=4
  if [ "$Nproc" -le "$Nsugg" ]; then let Nsugg=($Nproc+1)/2; fi
  echo "==================================================================="
  if [ "$Nproc" -eq "1" ]; then
    echo " Running the optimization process."
    echo " You can monitor the progress with"
    echo " $ tail -f proc*/output"
    echo " You can kill it with"
    echo " $ pkill -f main"
    echo " For more info, execute"
    echo " $ ./optimize.sh help"
  else
    echo " Running all $Nproc optimization processes at once in parallel."
    echo " You can monitor the progress with"
    echo " $ tail -f proc*/output"
    echo " You can kill all processes with"
    echo " $ pkill -f main"
    echo " You can restrict the number of parallel processes to, say $Nsugg, with"
    echo " $ ./optimize.sh Nparallel=$Nsugg"
    echo " For more info, execute"
    echo " $ ./optimize.sh help"
  fi
  echo "==================================================================="
fi

if [ $Nproc -lt $Ncpu ]; then
  Neven=$Nproc
else
  let Neven=($Nproc/$Ncpu)*$Ncpu
fi

function procQ {
  for ((i=$1; i<$Neven; i=i+$Ncpu)); do
    eval ${proc[$i]}
  done
  let a=$Neven+$1
  if [ $a -lt $Nproc ]; then
    eval ${proc[$a]}
  fi
}

for ((j=0; j<$Ncpu; j++)); do
  procQ $j &
done

#wait

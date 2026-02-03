#! /bin/bash

set -e

#Begin main list
Nmain=14
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
#End main list

for ((i=0; i<$Nmain; i++)); do
  echo ${array[$i]}
  eval ${array[$i]}
done


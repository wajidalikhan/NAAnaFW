#!/bin/bash

# Check if option is specified
if [[ "$1" == "" ]] 
then
    echo Select an input file
    exit
fi

inputfile=$1

python diffNuisances.py $inputfile  --vtol=0.001 -a -f text  | sed -e 's/!/ /g' -e 's/,/ /g' | sed -e 's/*/ /g' -e 's/,/ /g' | tail -n +2 > pulls.txt # -g pulls.root
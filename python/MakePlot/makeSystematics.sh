#!/usr/bin/env bash
declare -a arr=("jesUp" "jesDown" "btagDown" "btagUp" "puDown" "puUp" "lepDown" "lepUp" "jesUp" "jesDown" "mistagDown" "mistagUp")
for i in "${arr[@]}"
do
  echo "*** Producing histos for variable : $i *** "
  python makePlots.py -c muon -l 35.89 -r . --sys $i
done


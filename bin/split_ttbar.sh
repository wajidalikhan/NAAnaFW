#!/bin/sh
infile=TT.txt
path=work/files
rm -r $path/TT_*

lines=`awk 'END {print NR,"lines"}' $path/$infile`
echo 'Total No. of lines in file: '$infile : $lines 
#for i in $(seq 0 30 `awk 'END {print NR}' $path/$infile`) 
#  do
#  echo $i
  sed -n '1,   30  p' work/files/TT.txt > work/files/TT_1_30.txt
  sed -n '31,  60  p' work/files/TT.txt > work/files/TT_31_60.txt
  sed -n '61,  90  p' work/files/TT.txt > work/files/TT_61_90.txt
  sed -n '91,  120 p' work/files/TT.txt > work/files/TT_91_120.txt
  sed -n '121, 150 p' work/files/TT.txt > work/files/TT_121_150.txt
  sed -n '151, 180 p' work/files/TT.txt > work/files/TT_151_180.txt
  sed -n '181, 210 p' work/files/TT.txt > work/files/TT_181_210.txt
  sed -n '211, 240 p' work/files/TT.txt > work/files/TT_211_240.txt
  sed -n '241, 270 p' work/files/TT.txt > work/files/TT_241_270.txt
  sed -n '271, 300 p' work/files/TT.txt > work/files/TT_271_300.txt
  sed -n '301, 330 p' work/files/TT.txt > work/files/TT_301_330.txt
  sed -n '331, 360 p' work/files/TT.txt > work/files/TT_331_360.txt
  sed -n '361, 390 p' work/files/TT.txt > work/files/TT_361_390.txt
  sed -n '391, 420 p' work/files/TT.txt > work/files/TT_391_420.txt
  sed -n '421, 450 p' work/files/TT.txt > work/files/TT_421_450.txt
  sed -n '451, 480 p' work/files/TT.txt > work/files/TT_451_480.txt
  sed -n '481, 510 p' work/files/TT.txt > work/files/TT_481_510.txt
  sed -n '511, 540 p' work/files/TT.txt > work/files/TT_511_540.txt
  sed -n '541, 570 p' work/files/TT.txt > work/files/TT_541_570.txt
  sed -n '571, 600 p' work/files/TT.txt > work/files/TT_571_600.txt
  sed -n '601, 630 p' work/files/TT.txt > work/files/TT_601_630.txt
  sed -n '631, 660 p' work/files/TT.txt > work/files/TT_631_660.txt
  sed -n '661, 690 p' work/files/TT.txt > work/files/TT_661_690.txt
  sed -n '691, 720 p' work/files/TT.txt > work/files/TT_691_720.txt
  sed -n '721, 750 p' work/files/TT.txt > work/files/TT_721_750.txt
  sed -n '751, 780 p' work/files/TT.txt > work/files/TT_751_780.txt
  sed -n '781, 810 p' work/files/TT.txt > work/files/TT_781_810.txt
  sed -n '811, 840 p' work/files/TT.txt > work/files/TT_811_840.txt
  sed -n '841, 870 p' work/files/TT.txt > work/files/TT_841_870.txt
  sed -n '871, 900 p' work/files/TT.txt > work/files/TT_871_900.txt
  sed -n '901, 930 p' work/files/TT.txt > work/files/TT_901_930.txt
  sed -n '931, 960 p' work/files/TT.txt > work/files/TT_931_960.txt
  sed -n '961, 1000 p' work/files/TT.txt > work/files/TT_961_1000.txt
#  done
SingleTopAnalysis TT_1_30 work/files/TT_1_30.txt fullhadronic cat2 noSys noSync MC
SingleTopAnalysis TT_31_60 work/files/TT_31_60.txt fullhadronic cat2 noSys noSync MC

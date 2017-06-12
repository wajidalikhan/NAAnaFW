#!/bin/sh
#rm -r list.txt
#echo "Operation expired -- or -- No servers are available to read the file"
#echo "Job in Getting the list of jobId's from CERN Batch Cluster ..."

#path=/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/work/files
path=$CMSSW_BASE/Analysis/NAAnaFW/bin/work/files
#bjobs -r -o JOB_NAME | grep -v JOB_NAME | cut -d ':' -f2 | awk -F "_muon" '{print $1}' > list.txt
bjobs -r -o JOB_NAME | grep -v JOB_NAME | cut -d ':' -f2 | awk -F "_muon" '{print $1}' | awk '{print $0 ".txt"}' > list.txt

for i in `cat list.txt`;
  do
    if [[ -f $path/$i ]] ; then
          echo "The file $path/$i exists."
          #sed -i -e 's/xrootd-cms.infn.it/xrootd.ba.infn.it/g' /$path/$i
          #sed -i -e 's/xrootd.ba.infn.it/xrootd.ba.infn.it/g' /$path/$i
    else
          echo "The file $path/$i does not exist."
    fi
  done


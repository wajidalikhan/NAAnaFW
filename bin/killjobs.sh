#!/bin/sh
rm -r listjobs.txt
echo "Getting the list of jobId's from CERN Batch Cluster ..."
bjobs -rp -o jobid | grep -v ^JOBID | cut -d ':' -f2 > listjobs.txt
#bjobs -r -o jobid | grep -v ^JOBID | cut -d ':' -f2 > listjobs.txt

for i in `cat listjobs.txt`;
    do
      bkill $i;
    done


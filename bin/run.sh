#!/bin/bash

nbeg=$1
nend=$2

for ((j=${nbeg};j<${nend};j++))
do
    #!/bin/sh                                                                                         
    sbatch slurmJob-smarts.pbs ${j}
    sleep 0.5
done

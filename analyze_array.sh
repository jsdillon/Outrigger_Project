#! /bin/bash
#$ -V
#$ -S /bin/bash

cd /nfs/eor-00/h1/jsdillon/Outrigger_Project/
which python
source activate ENV
which python
echo $LD_LIBRARY_PATH

python Outrigger_Mapmaking.py $((SGE_TASK_ID - 1)) > ./Results/Logs/log_$((SGE_TASK_ID - 1)).txt

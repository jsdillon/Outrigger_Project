#! /bin/bash
#$ -V
#$ -S /bin/bash

exec bash
cd /nfs/eor-00/h1/jsdillon/Outrigger_Project/
which python
source activate ENV
which python

python Outrigger_Mapmaking.py $1 > ./Results/Logs/log_$1.txt

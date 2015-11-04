#!/bin/sh

echo "Process start number?"
read nStart

echo "Process end number?"
read nEnd

qsub -t $((nStart + 1))-$((nEnd + 1)) -V -l h_rt=48:00:00,h=eor-04,h_vmem=3G -pe chost 1 ./analyze_array.sh

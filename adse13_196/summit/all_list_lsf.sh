#!/bin/bash
#BSUB -P CHM137
#BSUB -W 02:00
#BSUB -nnodes 1
#BSUB -o job%J.out
#BSUB -e job%J.err
cd $WORK

echo "jobstart $(date)";pwd;ls
./all_list.sh
echo "jobend $(date)";pwd

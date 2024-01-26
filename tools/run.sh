#!/bin/bash
#SBATCH -p moh1
#SBATCH --mem=40G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00-15:00:00
#SBATCH -o run.sh_%j.out
#SBATCH -e run.sh_%j.err

tic=`date "+%Y-%m-%d %H:%M:%S"`;
printf "start at $tic\n\n";

make;
./sparse -I data/library -m 1 -q data/query;

toc=`date "+%Y-%m-%d %H:%M:%S"`;
printf "finish at $toc\n\n";

#!/bin/bash

echo "=================================== st.sh ==================================="  

# rm  ~/_scratch/log* 
rm ~/_scratch/slurm* ~/_scratch/main
scancel -u asashabokov98_607
mpicc ~/main.c -o ~/_scratch/main
sbatch -p test -n 16 ompi ~/_scratch/main ~/_scratch/input.txt

while [ ! -f ~/_scratch/slurm* ]
do
	sleep 1
done

# echo "Exists"
while [ ! -s ~/_scratch/slurm* ]
do
	sleep 1
done

cat ~/_scratch/slurm* | less 

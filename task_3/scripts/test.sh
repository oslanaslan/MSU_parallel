#!/bin/bash

for ((i=1; i < 33; i++))
do 
	for ((j=0; j < 10; j++))
	do 
		rm ~/_scratch/slurm*
		p=$(( $i * 8 )) 
		echo "$p  $j"
		cd ~/_scratch/
		sbatch -p test -n $p ompi ~/_scratch/main ~/_scratch/input_$j.txt
		while [ ! -f ~/_scratch/slurm* ]
		do 
			sleep 1
		done 
		while [ ! -s ~/_scratch/slurm* ]
		do
			sleep 2
		done
	done
done


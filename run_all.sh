#!/bin/bash

#for s in 1; do
#	for u in 0.00 0.02 0.04 0.06; do 
#		for j in 0.00 0.02 0.04 0.06; do 
#			echo ./SelfishDNA -o R5_Spatial_ABpulse_Sweep_j${j}_u${u}_s${s} -Size 150 -Scale 1 -Diff 2 -Jumprate ${j} -Uptakerate ${u} -ABinf 0.025 -ABstart 100000 -ABinterval 25000 -ABduration 12500 -Seed ${s} -x; 
#		done; 
#	done;
#done; 

size=100
for s in 1; do
	for b in 0.0 0.5 1.0; do
		for u in 0.01 0.02 0.03; do 
			for j in 0.03; do	
				echo ./SelfishDNA -o R8_Spatial_ABconst_Sweep_j${j}_u${u}_s${s}_b${b} -Size $size -Scale 1 -Diff 2 -Jumprate ${j} -Uptakerate ${u} -ABinf 0.0125 -ABstart 200000 -ABinterval 50000 -ABduration 99915000 -Seed ${s} -x -BreakChance ${b}; 
				echo ./SelfishDNA -o R8_Mixed_ABconst_Sweep_j${j}_u${u}_s${s}_b${b} -Size $size -Scale 1 -Diff 2 -Jumprate ${j} -Uptakerate ${u} -ABinf 0.0125 -ABstart 200000 -ABinterval 50000 -ABduration 99915000 -Seed ${s} -x -MixPop -BreakChance ${b}; 
			done; 
		done;
	done;
done;



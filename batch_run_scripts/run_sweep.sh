
size=150
t=100000
cost=0.005
for s in {11..15}; do
#	for phi in $(seq 0.1 0.1 1.0); do
	for phi in 0.85; do
		#for n in 0 1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80; do
		for n in 10 80; do
       			#echo ./SelfishDNA -o Sweep_MixPop_Phi${phi}_Non${n}_s${s} -Size $size -Diff 0.02 -Phi $phi -Jumprate 0.01 -Uptakerate 0.01 -Seed $s -BreakChance 1.00 -Cost $cost -Scale 2 -HK 10 -Noncoding $n -MaxTime $t -MixPop -x;
	       		#echo ./SelfishDNA -o Sweep_Space_Phi${phi}_Non${n}_s${s} -Size $size -Diff 0.02 -Phi $phi -Jumprate 0.01 -Uptakerate 0.01 -Seed $s -BreakChance 1.00 -Cost $cost -Scale 2 -HK 10 -Noncoding $n -MaxTime $t -x;
       			#echo ./SelfishDNA -o Sweep_MixDNA_Phi${phi}_Non${n} -Size $size -Diff 0.02 -Phi $phi -Jumprate 0.01 -Uptakerate 0.01 -Seed $s -BreakChance 0.00 -Cost $cost -Scale 2 -HK 10 -Noncoding $n -MaxTime $t -MixDNA -x;
       			echo ./SelfishDNA -o Log_ext_200_notims_Non${n}_s${s} -Size $size -Diff 0.01 -Phi $phi -Jumprate 0.02 -Uptakerate 0.01 -Seed $s -BreakChance 0.00 -Cost $cost -Scale 2 -HK 10 -Noncoding $n -MaxTime $t ;
		done;
	done;
done;


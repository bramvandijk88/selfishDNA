
for s in {1..10}; do 
	#echo ./SelfishDNA -o Log_Space_200_Anctrace_notimshgt_s${s} -Size 100 -Diff 0.02 -Phi 0.85 -Jumprate 0.02 -Uptakerate 0.01 -Seed $s -BreakChance 0.00 -Cost 0.005 -Scale 1 -HK 10 -Noncoding 40 -MaxTime 100000;
        #echo ./SelfishDNA -o Log_Space_200_Anctrace_nohgt_s${s} -Size 100 -Diff 0.02 -Phi 0.85 -Jumprate 0.02 -Uptakerate 0.00 -Seed $s -BreakChance 1.00 -Cost 0.005 -Scale 1 -HK 10 -Noncoding 40 -MaxTime 100000 -noHGT;
        echo ./SelfishDNA -o Log_nohgtlater_200_timshgt_s${s} -Size 150 -Diff 0.02 -Phi 0.85 -Jumprate 0.02 -Uptakerate 0.01 -Seed $s -BreakChance 0.50 -Cost 0.005 -Scale 1 -HK 10 -Noncoding 30 -MaxTime 200000;
        #echo ./SelfishDNA -o Log_Space_200_Anctrace_sex_s${s} -Size 100 -Diff 0.02 -Phi 0.85 -Jumprate 0.02 -Uptakerate 0.00 -Seed $s -BreakChance 1.00 -Cost 0.005 -Scale 1 -HK 10 -Noncoding 40 -MaxTime 100000 -noHGT -Sexual -Non 50;
        #echo ./SelfishDNA -o Log_Space_200_Anctrace_timsnocost_s${s} -Size 100 -Diff 0.02 -Phi 0.85 -Jumprate 0.02 -Uptakerate 0.01 -Seed $s -BreakChance 1.00 -Cost 0.000 -Scale 1 -HK 10 -Noncoding 40 -MaxTime 100000 -noHGT -Sexual;
done;


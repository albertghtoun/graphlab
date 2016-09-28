#!/bin/bash
out="outfile.txt"

echo -n "gather%%" > $out
for file in `ls input/*.txt`
do
	echo -n " $file" >> $out
done
echo >> $out

for gp in `seq 10`
do
	echo `echo "0.1 * $gp" | bc` 
	echo -n `echo "0.1 * $gp" | bc` >> $out
	for file in `ls input/*.txt`
	do
		echo "    $file"
		./debug/apps/pagerank/pagerank --infile $file --gather_percentage `echo "0.1 * $gp" | bc` 2>/dev/null | grep "runtime:" | awk '{printf(" %.2f"), $4}' >> $out
	done
	echo >> $out
done

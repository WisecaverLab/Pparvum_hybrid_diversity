#!/bin/bash
#cd ./02_orthofinder_oneoutgroup/Results_Mar07/Single_Copy_Orthologue_Sequences/


files=$( ls hybrid_orthofinder/17_orthofinder/Results_Jan03/Single_Copy_Orthologue_Sequences )
cd hybrid_orthofinder/scripts/scog_100_lists

#for file in $files

i=0
for file in $files
do
	filenum=$(( $i/100 ))
	echo $file >> scoglist_${filenum}.txt
	i=$(( $i + 1 ))
done

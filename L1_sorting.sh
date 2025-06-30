#!/bin/bash

while read L1; do
        id=$(echo $L1 | cut -f4 -d' ')
	awk '{if ($4 > 2) print $0 "\t" FILENAME}' */*$id.bed | sed 's/\//\t/g' | cut -f1,2,3,4,6 | sed 's/.bed//g' | sort -k1,1 -k2,2n > L1_based/$id"_all.bed"
	bedtools window -w 2000 -v -b <(grep $id near_reagion_with_id.bed | awk '{if ($4 < 3) print $0}') -a <( bedtools merge -d 2000 -i L1_based/$id"_all.bed" -o collapse,count -c 5) > L1_based/$id"_merged.bed"
done <$1

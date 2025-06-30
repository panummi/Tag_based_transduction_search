#!/bin/bash

#  Created on Tuesday, 21 December  2021

# For debugging
#set -o verbose

# Die on unset variables
#set -o nounset
# Die on errors
#set -o errexit
# Die if any part of a pipe fails
#set -o pipefail


cram="$1"
test -e "${cram}"

MAX_THREADS=5

seqs=$2

TEMPDIR=$(mktemp --directory)
        function _cleanup {
                echo "Cleaning up." >&2
        	ls -ltr "${TEMPDIR}" >&2
        	rm -r "$TEMPDIR"
        }
trap _cleanup EXIT

pattern=""

##Iterate over all sequences and add them and the reverse complement to var pattern

while read te; do
	seq1=$(echo $te | cut -f2 -d' ')
	id=$(echo $te | cut -f1 -d' ')
	echo $id
	seq2=$(echo $seq1 | tr ACGTacgt TGCAtgca | rev)


	pattern=$pattern"|${seq1}|${seq2}"

done <$seqs

pattern=$(echo $pattern | sed 's/|//1')

echo "Starting to grep ${cram} at $(date) on $(uname --all)" >&2

# Using ripgrep search for all the different sequences to search the bam files
GREP="rg"

/usr/bin/time -v samtools view -@ ${MAX_THREADS} "$cram" | ${GREP} "${pattern}" >${TEMPDIR}/samish
echo "Done grepping at $(date)" >&2

mkdir -p $(basename $cram .cram)


while read te; do
#        echo $te
        seq1=$(echo $te | cut -f2 -d' ')
        id=$(echo $te | cut -f1 -d' ')
	echo $id
        seq2=$(echo $seq1 | tr ACGTacgt TGCAtgca | rev)
        outbam=$(basename $cram .cram)"/"$(basename $cram .cram).$id.bam
        outbed=$(basename $cram .cram)"/"$(basename $cram .cram).$id.bed

	#search for one sequence at the time from a temporary file where all matches from all the used sequences are

	pattern="${seq1}|${seq2}"
	${GREP} "${pattern}" ${TEMPDIR}/samish > ${TEMPDIR}/${id}_samish || echo Zero matches found

	#make a bed file to use in search for the anchor reads

	awk -v OFS=$'\t' '{print $3,$4-1,$4+2;print ($7=="="?$3:$7),$8-1,$8+1;}' ${TEMPDIR}/${id}_samish | awk '{if ($1 != "*") print $0}' |
    	  LC_ALL=C sort --key=1,1 --key=2,2n |
    	  bedtools merge -d 10 -i - >${TEMPDIR}/${id}_read_regions.bed


	cut -f1 ${TEMPDIR}/${id}_samish | sort -u > ${TEMPDIR}/${id}_read_names

	#extract the read pairs, anchor reads, from original cram file using read names and a bed file

	/usr/bin/time -v samtools view -M -@ ${MAX_THREADS} -L ${TEMPDIR}/${id}_read_regions.bed -o "${outbam}" --qname-file ${TEMPDIR}/${id}_read_names "${cram}"

	#merge reads based on the anchor loci to create insertion calls

	n=$(bedtools bamtobed -i "${outbam}" | awk '$5>36' | wc -l)
	if (( n > 0 )); then
		/usr/bin/time -v bedtools bamtobed -i "${outbam}" | awk '$5>36' | bedtools merge -d 1000 -c 4 -o count >"${outbed}"
	fi
	echo "All done at $(date)" >&2

done <$seqs

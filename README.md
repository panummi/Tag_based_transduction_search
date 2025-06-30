# Tag_based_transduction_search
Detecting short sequences originating downstream of source L1 elements and identifying 3'transducions in short-read Illumina data
Usage:

count_unique.py

Used to count the number of times a tag-sequence or sequences like it are present in the used reference genome.

use: python count_unique.py -i input_of_sequences.fa -o output_name -r reference_genome.fa -b bed_file_for_output.bed

tag_based_search.sh

Used to find tag sequences in discordant reads:

use: bash tag_based_search.sh alignment_file.cram tag_sequences.csv 

L1_sorting.sh 

Used to merge the calls by source L1s:

use: bash source L1s tag_sequences.csv 

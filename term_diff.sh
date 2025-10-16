#!/bin/bash
# To run :
# chmod +x term_diff.sh
# ./term_diff.sh
GENOME_FILE="/path/to/complete/genome.fasta"
diff <(seqkit subseq -r 1:1000 ${GENOME_FILE} | grep -v ">") \
     <(seqkit subseq -r -1000:-1 ${GENOME_FILE} | grep -v ">")

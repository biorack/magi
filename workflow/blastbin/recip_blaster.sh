#!/bin/bash

#to be put in the same folder as the blast binaries

#number of cpus
CPU=$1
#blast database path
DB=$2
#working dir
CWD=$3
#repo location
REPO=$4

i=$CPU

BLASTP="$REPO/workflow/blastbin/blastp"

while [ $i -gt 0 ]
do
	i=`expr $i - 1`
	SEQ="$CWD/tmp_seq_$i.faa"
	OUT="$CWD/tmp_out_blasted_$i.txt"
	START=$(date +%s.%N)
	$BLASTP -query $SEQ -db $DB -outfmt '10 qacc sacc qcovs length ppos evalue bitscore' -evalue 1 -max_target_seqs 10 > $OUT &
done

wait

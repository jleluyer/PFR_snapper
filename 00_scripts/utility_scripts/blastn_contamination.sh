#!/bin/bash
#PBS -N blastn
#PBS -o blastn_valid.out
#PBS -l walltime=20:00:00
#PBS -l mem=30g
#PBS -l ncpus=16
#PBS -q omp
#PBS -r n



cd $PBS_O_WORKDIR

# baslt variables
OUTFMT="\"6 qseqid qlen sseqid sallacc slen qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore\""

# Global variables
TRANSCRIPTOME="transcriptome.fa"



bacteria="refseq/refseq.batceria.fna"

#makeblastdb -in $bacteria -input_type 'fasta' -dbtype 'nucl'

#cat $TRANSCRIPTOME | parallel -j 16 -k --block 10k --recstart '>' --pipe blastn -db $bacteria -query - -outfmt $OUTFMT -max_target_seqs 1 -evalue 1e-10 >blastn_bacteria_e-4.txt



viral="refseq/refseq.viral.fna"

#makeblastdb -in $viral -input_type 'fasta' -dbtype 'nucl'

#cat $TRANSCRIPTOME | parallel -j 16 -k --block 10k --recstart '>' --pipe blastn -db $viral -query - -outfmt $OUTFMT -max_target_seqs 1 -evalue 1e-10 >blastn_virus_e-10.txt

fungi="refseq/refseq.fungi.rna.fna"

#makeblastdb -in $fungi -input_type 'fasta' -dbtype 'nucl'

#cat $TRANSCRIPTOME | parallel -j 16 -k --block 10k --recstart '>' --pipe blastn -db $fungi -query - -outfmt $OUTFMT -max_target_seqs 1 -evalue 1e-10 >blastn_fungi_e-10.txt


plasmid="refseq/refseq.plasmid.rna.fna"

#makeblastdb -in $plasmid -input_type 'fasta' -dbtype 'nucl'

cat $TRANSCRIPTOME | parallel -j 16 -k --block 10k --recstart '>' --pipe blastn -db $plasmid -query - -outfmt $OUTFMT -max_target_seqs 1 -evalue 1e-10 >blastn_plasmid_e-10.txt

archaea="refseq/refseq.archaea.fna"
#makeblastdb -in $archaea -input_type 'fasta' -dbtype 'nucl'

cat $TRANSCRIPTOME | parallel -j 16 -k --block 10k --recstart '>' --pipe blastn -db $archaea -query - -outfmt $OUTFMT -max_target_seqs 1 -evalue 1e-10 >blastn_archaea_e-10.txt


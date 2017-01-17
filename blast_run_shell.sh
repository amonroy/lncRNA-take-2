#!/bin/bash

#SBATCH --job-name=blastn
#SBATCH --ntasks=1
#SBATCH --time=7-1:00:00
#SBATCH --mem=30000
#SBATCH --cpus-per-task=12
#SBATCH -o blastn_%A.out
#SBATCH -e blastn_%A.err

module add blast

for i in /proj/cdjones_lab/anais_proj/blast2/*.fasta
do
    ##echo "This is me checking:"
    ##echo $i
    name=$(basename $i .fasta)
    ###eval $name
    ##echo "this is me seeing if command worked"
    ##echo "$i -> $name"
    blastn -db $i -query /proj/cdjones_lab/anais_proj/dmel/dmel-all-transcript-r6.13.fasta  -out /proj/cdjones_lab/anais_proj/blast_out/$name+_x_dmel.blstn -outfmt 6 -num_threads 12



done


##do blastn -db -query -out -outfmt 6 -num_threads 12


##for i in /proj/cdjones_lab/anais_proj/blast2/*.fasta;
##do makeblastdb -in "$i"  -dbtype nucl -parse_seqids;
##done

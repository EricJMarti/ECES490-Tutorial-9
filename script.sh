#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -M ejm335@drexel.edu
#$ -P nsftuesPrj
#$ -pe shm 16-32
#$ -l h_rt=48:00:00
#$ -q all.q@@amdhosts

. /etc/profile.d/modules.sh
module load shared
module load proteus
module load sge/univa
module load gcc/4.8.1

## Convert SRA files to FASTA files
#PATH=/mnt/HA/groups/nsftuesGrp/.local/bin:$PATH

#fastq-dump --fasta *.sra

## Run FragGeneScan on FASTA files
datafiles=(SRR492065 SRR492066 SRR492182 SRR492183 SRR492184 SRR492185 SRR492186 SRR492187 SRR492188 SRR492189 SRR492190 SRR492191 SRR492192 SRR492193 SRR492194 SRR492195 SRR492196 SRR492197)

FRAGPATH=/home/ejm335/fraggenescan-code

for file in ${datafiles[@]}
do
$FRAGPATH/run_FragGeneScan.pl -genome=./$file.fasta -out=./output/$file -complete=0 -train=454_30
done

#!/bin/bash

#SBATCH --time=100:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=gatk_merge_edit

module purge

module load gatk/4.2.0.0
module load picard/2.23.8
module load bwa/intel/0.7.17
module load samtools/intel/1.14/
module load bamtools/intel/2.5.1

REF=/scratch/cb4097/Sequencing/Reference/*.fna

echo $REF

#echo $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}
#echo "HTG3TDMXY is the new merged file name"

#samtools merge HNGLVDRXY_gR_merged.bam gR*.bam
#samtools merge HNGLVDRXY_${1}_merged.bam $2

samtools merge AllCuSO4 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} ${21} ${22} ${23} ${24} ${25} ${26} ${27} ${28} ${29} ${30} ${31} ${32} ${33} ${34} ${35} ${36} ${37} ${38} ${39} ${40} ${41} ${42} ${43} ${44} ${45} ${46} ${47} ${48} ${49} ${50} ${51} ${52} ${53} ${54} ${55} ${56} ${57} ${58} ${59} ${60} ${61} ${62} ${63} ${64} ${65} ${66} ${67} ${68} ${69} ${70}

samtools index AllCuSO4

#split the files
bamtools split -in AllCuSO4 -reference

#Don't run this; do it in parallel because of timing:
#gatk HaplotypeCaller -I $2 -R $REF -ploidy 1 -O HNGLVDRXY_${1}_merged.vcf

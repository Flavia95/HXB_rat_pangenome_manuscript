#!/bin/bash
# Usage: small_evaluation.sh <SDF of reference> <chr of interest>

#SAMPLE=$1
thread=10
SDF=$1
chr=$2

rtg vcfeval -t $SDF -b truth_vcf/$2/deepvariant_$2_33_samples_reh.mRatBN7.2.gvcf.gz.max50.SHROlaIpcv.$2.vcf.gz -c graphs/$2.pan+ref/tmp.2.max50.SHROlaIpcv.$2.vcf.gz --squash-ploidy -e ucsc_easy_hard/easy.$2.bed -T ${thread} -o vcfeval_easy
rtg vcfeval -t $SDF -b truth_vcf/$2/deepvariant_$2_33_samples_reh.mRatBN7.2.gvcf.gz.max50.SHROlaIpcv.$2.vcf.gz -c graphs/$2.pan+ref/tmp.2.max50.SHROlaIpcv.$2.vcf.gz --squash-ploidy -e ucsc_easy_hard/hard_$2.bed -T ${thread} -o vcfeval_hard

# snps and indels only easy region
rtg vcfeval -t $SDF -b truth_vcf/$2/SHROlaIpcv.JT.snps.vcf.gz -c graphs/$2.pan+ref/SHROlaIpcv.pggb.snps.vcf.gz --squash-ploidy -e ucsc_easy_hard/easy.$2.bed -T ${thread} -o vcfeval_snps
rtg vcfeval -t $SDF -b truth_vcf/$2/SHROlaIpcv.JT.indels.vcf.gz -c graphs/$2.pan+ref/SHROlaIpcv.pggb.indels.vcf.gz -e ucsc_easy_hard/easy.$2.bed -T ${thread} --squash-ploidy  -o vcfeval_indels

# hard region
rtg vcfeval -t $SDF -btruth_vcf/$2/SHROlaIpcv.JT.snps.vcf.gz -c ucsc_easy_hard/hard_$2.bed -T ${thread} --squash-ploidy -o vcfeval_snps_hard
'''
#cd vcfeval_easy/SHROlaIpcv
#echo sample tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t' > statistics.easy.2000.F100.tsv
#grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' >> statistics.easy.2000.F100.tsv
#cd ..
#cd vcfeval_snps/SHROlaIpcv
#echo sample tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t' > statistics.snps.easy.tsv
#grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' >> statistics.snps.easy.tsv
#cd ..
#cd vcfeval_indels/SHROlaIpcv
#echo sample tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t' > statistics.indels.easy.tsv
#grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' >> statistics.indels.easy.tsv
#cd ..

#cd vcfeval_hard/SHROlaIpcv
#echo sample tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t' > statistics.hard.2000.F100.tsv
#grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' >> statistics.hard.2000.F100.tsv
#cd ..

#cd vcfeval_snps_hard/SHROlaIpcv
#echo sample tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t' > statistics.snps.hard.tsv
#grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' >> statistics.snps.hard.tsv
#cd ..
'''

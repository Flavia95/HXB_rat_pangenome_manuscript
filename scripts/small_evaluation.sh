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

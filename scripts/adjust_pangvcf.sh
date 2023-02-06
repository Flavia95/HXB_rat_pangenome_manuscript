#!/bin/bash
# Usage: adjust_pangvcf.sh <vcf file obtained by pggb>

input=$1

#1. Remove Ns sites from VCF                                                                                                        
awk '$5  $4 !~ /N/ {print$ALL}' $1 > pan+ref.fa.gz.smooth.fix.REF.noN.vcf                                                                                                                                                                                                     
#2. Adjust header and chr column                                                                                                                                     
bcftools annotate --rename-chrs chr_name_conv.txt pan+ref.fa.gz.smooth.fix.REF.noN.vcf | bgzip > tmp.2.vcf.gz
#3. Normalize,decompose,consider only small variants and split by snps and indels
for v in `ls tmp.2.vcf.gz`;
do
	for s in `cat rename.txt`;
	do
		bash vcf_preprocess.pang.sh ${v} ${s} 50
	done
done

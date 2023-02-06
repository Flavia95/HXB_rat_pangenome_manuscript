#!/bin/bash


# From https://github.com/wwliao/pangenome-utils/blob/main/preprocess_vcf.sh
# Usage: preprocess_vcf.sh <VCF file> <sample name> <max variant size>

VCF=$1
FNAME=$(basename $VCF)
PREFIX=$(dirname $VCF)/"${FNAME%.vcf.gz}"
SAMPLE=$2
CHROMS='chr1'
MAXSIZE=$3
REF="/lizardfs/flaviav/rat/script/rn7.fa"
MEM="10G"

#Filter VCFs
bcftools=/home/flaviav/.guix-profile/bin/bcftools
$bcftools view -a -s ${SAMPLE} -Ou ${VCF} \
	| $bcftools norm -f ${REF} -c s -m - -Ou \
	| $bcftools view -e 'GT="ref" | GT~"\."' -f 'PASS,.' -Ou \
	| $bcftools sort -m ${MEM} -T bcftools-sort.XXXXXX -Ou \
	| $bcftools norm -d exact -Oz -o ${PREFIX}.${SAMPLE}.norm.vcf.gz \
	&& $bcftools index -t ${PREFIX}.${SAMPLE}.norm.vcf.gz \
	&& $bcftools view -e "STRLEN(REF)>${MAXSIZE} | STRLEN(ALT)>${MAXSIZE}" \
	-r ${CHROMS} -Oz -o ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1.vcf.gz \
	${PREFIX}.${SAMPLE}.norm.vcf.gz \
	rm ${PREFIX}.${SAMPLE}.norm.vcf.gz*


#for F in *.vcf.gz ; do   tabix -f -p vcf ${F}  ; done

#Filter for SNPs and INDELs
$bcftools view -v snps ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1.vcf.gz -O z -o ${SAMPLE}.pggb.snps.vcf.gz
tabix ${SAMPLE}.pggb.snps.vcf.gz
$bcftools view -v indels ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1.vcf.gz -O z -o ${SAMPLE}.pggb.indels.vcf.gz
tabix ${SAMPLE}.pggb.dec.vcf.gz

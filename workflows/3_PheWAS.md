## PheWAS analyses

#### 1. Preprocessing

Merge all variant calling only by the pangenome and JC by sample after rtg evaluation.
This is an example:
```
bcftools merge */fp.vcf.gz -o onlysnps_pg_allsamples.vcf.gz
bcftools merge */fn.vcf.gz -o onlysnps_jc_allsamples.vcf.gz
```

#### 2. PheWAS analyses

For the 12 validated variants with the PCR I run the PheWAS pipeline.

- I extracted only the useful informations.  
```
bcftools query -f '%CHROM\t%POS\t[ %GT]\n' 12valiPCR.vcf.gz -o 12.onlyvg.nomulti.gt.txt
```

- [generategeno.R](https://github.com/Flavia95/HXB_pangenome/blob/main/scripts/generategeno.R), create ".geno" file from 12.onlyvg.nomulti.gt.txt, convert to a matrix file of genotypes.
- [generatepheno.R](https://github.com/Flavia95/HXB_pangenome/blob/main/scripts/generatepheno.R), consider single marker of interest each time and run the PHEWAS analyses with subset phenotype file dowloaded with the API from GeneNetwork database.

2 of 12 markers validated by PCR called by vg-only present a phenotypes associations.

##### 3. Explore these two markers and the phenotypes 
- Ensembl Genome Browser, selecting  mRatBN7.2 reference genome.
- Rat Genome Database [RGD](https://academic.oup.com/genetics/article-abstract/224/1/iyad042/7080056?redirectedFrom=fulltext&login=false)

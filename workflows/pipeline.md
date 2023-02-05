# 10x technology-pggb analysis pipeline

#### 1. Data and preprocessing

Create these folders
```
mkdir assemblies
mkdir alignments
mkdir parts
```
To change the sequence names according to PanSN-spec, we use fastix:
```
for f in *_megabubble.fasta.gz; do                    
        sample_name=$(echo $f | cut -f 1 -d '_');  
        fastix -p "${sample_name}#1#" <(zcat $f)> $f._changeid.fa                
        bgzip *_changeid.fa
        samtools faidx *_changeid.fa.gz  
`rename "s/.fasta.gz.//" *fa.gz`
```   
#### 2. Sequence divergence
The sequence divergence for each set of chromosomes is high because we are using inbreeding assemblies, so we are using -p 98.

#### 3. Sequence partitioning
Partition the assembly contigs by chromosome by mapping each assembly against the reference genome. We use wfmash for the mapping.

```
zcat * fasta.gz > rat.fasta.gz
zgrep '_changeid.fasta.gz$' | cut -f 1-2 -d . | sort -V | uniq >haps.list
ref=rn7_supernova_megabubble_changeid.fasta.gz
for hap in $(cat haps.list);
do
    in=assemblies/$hap
    out=alignments/$hap.vs.ref.paf
    wfmash -t 48 -m -N -p 90 -s 20000 $ref $in >$out
done
```
#### 4. Sequence partitioning
Partition the assembly contigs by chromosome by mapping each assembly against the reference genome. We use wfmash for the mapping.

```
zcat * fasta.gz > rat.fasta.gz
zgrep '_changeid.fasta.gz$' | cut -f 1-2 -d . | sort -V | uniq >haps.list
ref=rn7_supernova_megabubble_changeid.fasta.gz
for hap in $(cat haps.list);
do
    in=assemblies/$hap
    out=alignments/$hap.vs.ref.paf
    wfmash -t 48 -m -N -p 90 -s 20000 $ref $in >$out
done
```
#### 5. Subset by chromosome:

```
(seq 20; echo X; echo Y) | while read i; do awk '$6 ~ "chr'$i'$"' $(ls alignments/*.vs.ref.paf | sort -V) | cut -f 1 | sort -V > parts/chr$i.contigs; done 
(seq 20;echo X;echo Y) | while read i; do xargs samtools faidx rat.fasta.gz < parts/chr$i.contigs > parts/chr$i.pan.fa; done 
(seq 20; echo X; echo Y; echo M) | while read i; do samtools faidx rn7_supernova_megabubble_changeid.fasta.gz rn7#1#chr$i > rn7_changeid.chr$i.fa && cat rn7_changeid.chr$i.fa chr$i.pan.fa > chr$i.pan+ref.fa && bgzip chr$i.pan+ref.fa && samtools faidx chr$i.pan+ref.fa.gz; done
```

This results in chromosome-specific FASTAs in `parts/chr*.pan+ref.fa.gz`

#### 6. Graph generation (all chromosomes)

We apply [pggb](https://github.com/pangenome/pggb) and vg for the variant calling
```
(seq 20;echo X;echo Y| tr ' ' '\n') | while read i; do sbatch -p workers -c 48 --wrap 'cd /scratch && /lizardfs/flaviav/rat/rat_paper/pggb_c1c3a1565604fc41f880bccd9f46d0a709f3e774 -i /lizardfs/flaviav/rat/pggb_paper/parts/chr'$i'.pan+ref.fa.gz -p 98 -s 2000 -n 32 -F 0.001 -k 79 -P asm5 -O 0.03 -G 4001,4507 -V rn7:# -t 48 -T 48 -o chr'$i'.pan+ref ; mv /scratch/chr'$i'.pan+ref '$(pwd); done
```
#### 7. Adjust "truth set" (JT)
`bcftools reheader --samples rename.txt deepvariant_chr$1_33_samples_mRatBN7.2.gvcf.gz -o deepvariant_chr$1_33_samples_reh.mRatBN7.2.gvcf.gz`

```
for v in `ls deepvariant_chr$1_33_samples_reh.mRatBN7.2.gvcf.gz`;
do
    for s in `cat rename.txt`;
    do
        bash  [vcf_preprocess.sh](script/vcf_preprocess.sh) ${v} ${s} 50
    done
done
```
#### 8. Adjust pangenomic set (VG)

[adjust_pangvcf.sh](script/adjust_pangvcf.sh) fileobtainedbypggb.vcf

#### 9. Graph evaluation
```
mkdir evaluation
```
Prepare the reference genome:
```
rtg format -o ref.sdf rn7.fa
```
Evaluation between the "truth set" and pangenomic set.

[small_evaluation.sh](script/small_evaluation.sh) ref.sdf chrofinterest

### Results

- Graph statistics using odgi stats:
- Number of small variants using vt peek and bcftools stats:
- Evaluation graphs using rtg:

We focalized our work only using chr12, we selected some variants called only by the pangenome and we validated these with the PCR.

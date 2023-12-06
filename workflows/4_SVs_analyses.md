## Structural Variants analyses

To adopt a pangenomic perspective and capture SVs that might be missed with traditional reference-based methods, 

1.  We applied several SVs calling pipelines. 
-  As a graph-based method, we used PGGB and vg to generate a comprehensive VCF file. Variants were normalized using [vcfwawe](https://github.com/vcflib/vcflib) and then filtered by keeping only variants >50 bp.
- As assembly-based methods, we applied 3 pipelines (pav, hall, svim) as described in [A draft human pangenome reference](https://www.nature.com/articles/s41586-023-05896-x). 

2.  To evaluate the concordance between all the approaches, we conducted an overlap analysis using [truvari](https://github.com/ACEnglish/truvari) (v.4.0) 
 ```
    $SURVIVOR merge $SAMPLE.sample_files 1000 1 1 1 0 50 $SAMPLE.merged.tmp.vcf
   ```
We merge using a maximum allowed distance of 1000 bp, as measured pairwise between breakpoints (begin1 vs begin2, end1 vs end2).
We only to report calls supported by 2 callers and they have to agree on the type (1) and on the strand (1) of the SV.
We only compare SVs that are major of 50bp long.
We removed variants with more than 10% of Ns in REF/ALT. 

3. To validate the SVs identified, we employed Adaptive Sampling, a nanopore sequencing-specific software-controlled enrichment method.

4. For the SV’s pangenome visualization we used [odgi](https://github.com/pangenome/odgi) and [Bandage]()
-  We extracted only the SVs in which we are interesting used the coordinates of the reference with ```odgi extract -i graph.og -r reference:start-end -o extracted.og```.
-   We visulized SVs extracted with ```odgi viz``` and ```Bandage```.

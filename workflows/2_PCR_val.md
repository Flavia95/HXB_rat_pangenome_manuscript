## PCR validation pipeline

Variants only found with the pangenome approach were randomly selected from Chr 12 for the SHR/OlaIpcv sample. 
Sequencing results were mapped against the target sequence using the blastn and mutations were manually examined for confirmation.

#### 1. Reformat Sanger File
[reformat_sanger.pl](https://github.com/Flavia95/HXB_pangenome/blob/main/scripts/reformat_sanger.pl)
 
#### 2. Merge files
[merge_metadata.r](https://github.com/Flavia95/HXB_pangenome/blob/main/scripts/rat_pangenome_chr12_snp_validation_merge_metadata.r)

#### 3. Add search informations 
[add_search_info.sh](https://github.com/Flavia95/HXB_pangenome/blob/main/scripts/addsearch_info.sh)

#### 4. Run blast
[generate_blastresults.sh](https://github.com/Flavia95/HXB_pangenome/blob/main/scripts/generate_blastresults.sh)

#### 5. Manually curation for confirmation

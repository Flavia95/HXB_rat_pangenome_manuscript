#!/bin/bash
tr -d '"' <rat_pangenome_chr12_snp_validation_all_data_search_2023-02-03.csv >rat_pangenome_chr12_snp_validation_all_data_search_edit_2023-02-03.csv #remove quotes in a csv

awk -F ',' '{print ">"$1"\n"$20 > $1".seq.fa"; print ">"$2"\n"$19 > $2".targetSeq.fa"}' rat_pangenome_chr12_snp_validation_all_data_search_edit_2023-02-03.csv

while read line
do
search=$(echo $line | awk -F ',' '{print $21}')
alt=$(echo $line | awk -F ',' '{print $16}')
ref=$(echo $line | awk -F ',' '{print $15}')
query=$(echo $line | awk -F ',' '{print $1}')
target=$(echo $line | awk -F ',' '{print $2}')
seq_file=$query".seq.fa"
target_file=$target".targetSeq.fa"
if [ -f "$target_file" ]; then
blastn -subject $target_file -query $seq_file > $query"."$target".blast.txt"
echo "Search: "$search >> $query"."$target".blast.txt"
echo "ALT: "$alt >> $query"."$target".blast.txt"
echo "REF: "$ref >> $query"."$target".blast.txt"
else
echo "Error: target file $target_file not found for $seq_file"
fi
done < rat_pangenome_chr12_snp_validation_all_data_search_edit_2023-02-03.csv

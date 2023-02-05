#!/bin/bash
sed 's/\ /,/g'  rat_pangenome_chr12_snp_validation_all_data_2023-02-03.tab > output.csv  #transform tab file into a csv file
start=145
end=149
while read line; do
  arr=($(echo $line | tr ',' ' '))
  search=${arr[19]:start-1:end-start+1}${arr[15]:1:1}
  echo $line,$search
done < output.csv > rat_pangenome_chr12_snp_validation_all_data_search_2023-02-03.csv

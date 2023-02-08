#!/bin/bash
#Transform tab file into a csv file and remove quotes
sed 's/\ /,/g' rat_pangenome_chr12_snp_validation_all_data_2023-02-03.tab | tr -d '"' > rat_pangenome_chr12_snp_validation_all_data_2023-02-03.edit.csv

# Read the input file line by line
while IFS=, read -r ID name Sequence Scale Purification Tm begin tag direction CHROM start end newID POS REF ALT QUAL ab1name seq targetSeq
do
  if [ "$direction" == "R" ]; then
    # Use the tr and rev commands to perform the substitution and reversal on targetSeq (generate the reverse complementary)
    updated_targetSeq=$(echo "$targetSeq" | tr ACGT TGCA | rev)

    # Output the updated line to a new file
    echo "$ID,$name,$Sequence,$Scale,$Purification,$Tm,$begin,$tag,$direction,$CHROM,$start,end,$newID,$POS,$REF,$ALT,$QUAL,$ab1name,$seq,$updated_targetSeq" >> output.csv
  else
    # Output the unchanged line to the output file if direction is not R (it is F)
    echo "$ID,$name,$Sequence,$Scale,$Purification,$Tm,$begin,$tag,$direction,$CHROM,$start,end,$newID,$POS,$REF,$ALT,$QUAL,$ab1name,$seq,$targetSeq" >> output.csv
  fi
done < rat_pangenome_chr12_snp_validation_all_data_2023-02-03.edit.csv

#Removes the first row of the CSV file and adds a header row to the CSV file,
tail -n +2 output.csv > tmp.csv

#Performs a search operation on each line of the CSV file in the refseq (targetseq) and adds the result to the end of each line
echo "ID,name,Sequence,Scale,Purification,Tm,begin,tag,direction,CHROM,start,end,newID,POS,REF,ALT,QUAL,ab1name,seq,targetSeq,search" > rat_pangenome_chr12_snp_validation_all_data_search_2023-02-03.csvstart=146
while read line; do
arr=($(echo $line | tr ',' ' '))
if [ "${arr[8]}" == "F" ]; then
end=149
else
end=150
fi
search=${arr[19]:start-1:end-start+1}${arr[15]:1:1}
echo $line,$search >> rat_pangenome_chr12_snp_validation_all_data_search_2023-02-03.csv
done < tmp.csv

rm tmp.csv

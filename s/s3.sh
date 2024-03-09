#!/usr/bin/bash


configPath=$1
. $configPath


#outer functions and data paths:
make_genotype_pickle=${postERVcaller}/s/make_genotype_pickle.py
make_filtered_dataset_vcfRecords=${postERVcaller}/s/make_filtered_dataset_vcfRecords.py


#shell argument in
dataset=fltSetO1
infile_feature=${cohort}_ERVcaller
infile_feature2=${infile_feature}_${dataset}_step1

#out basenames
outfile_feature1=${infile_feature}_${dataset}
outfile_feature2=${outfile_feature1}_represenPos
outfile_feature3=${outfile_feature1}_represenID
outfile_feature4=${outfile_feature3}_clustermergedReGT

#other in files
in_table1=${infile_feature2}_pos_cluster_info.txt
in_table2=${infile_feature2}_clustermerged_cluster_info.txt
in_addIDvcf=${infile_feature}_addID.vcf
in_list_represenIDtoDel=${outfile_feature3}_toDelete_list.txt
in_clustermergedVcf=${infile_feature2}_represenID_clustermergedReGT.vcf
in_clustermergedVcfRec=${in_clustermergedVcf}_records


#For-PopGen out files
outfile4=${outfile_feature2}_pos_geno.pkl
f_clustermergedVcf=${outfile_feature4}.vcf
f_clustermergedVcfRec=${f_clustermergedVcf}_records
f_pos_geno=${outfile_feature4}_pos_geno.pkl


#other intermediate/out files
f_list=${outfile_feature2}_list.txt
f_list2=${outfile_feature3}_list.txt


#---make genotype matrix bg
#get subset records
#split header lines and variant records
f_header=${in_clustermergedVcf}_header
f_records=${in_clustermergedVcfRec}
awk -v f_header=$f_header -v f_records=$f_records '{if ($0 ~ /^##/){print $0 > f_header} else {print $0 > f_records}}' ${in_clustermergedVcf}\

#final filtering
${python3} ${make_filtered_dataset_vcfRecords} ${in_clustermergedVcfRec} ${in_table1} ${in_table2} ${f_clustermergedVcfRec} \
${in_list_represenIDtoDel} ${f_list} ${f_list2} || \
(echo "Failed: make_filtered_dataset_vcfRecords" && exit)
cat $f_header ${f_clustermergedVcfRec} > ${f_clustermergedVcf}
#---make genotype matrix ed
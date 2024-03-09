#!/usr/bin/bash


configPath=$1
. $configPath




#shell argument in
dataset=preflt
infile_feature=${cohort}_ERVcaller

#out basenames
outfile_feature1=${infile_feature}_${dataset}

#other in files
in_table3=${infile_feature}_preflt_posID_list.txt

#out files
f_addIDvcf=${infile_feature}_addID.vcf


#---generate some useful intermediate files bg
#add position ID to the original ERVcaller out merge vcf for later subsetting the vcf file.
readd=0
if ! [ -e $f_addIDvcf ]; then
	readd=1
else
	read -a alen1 <<< $(wc -l ${infile_feature}_merge.vcf);len1=${alen1[0]}
	read -a alen2 <<< $(wc -l $f_addIDvcf);len2=${alen2[0]}
	readd=$((${len1}-${len2}))
fi
if ! [ $readd -eq 0 ]; then
	awk '{if ($0 !~ /^#/) {new5=$5; sub(/<INS_MEI:/,"",new5); sub(/>/,"",new5); new3=$1":"$2":"new5;\
	sub(/\./,new3,$0)}; print $0}' ${infile_feature}_merge.vcf > $f_addIDvcf
fi

#split header lines and variant records
f_header=${f_addIDvcf}_header
f_records=${f_addIDvcf}_records
if ! [ -e $f_records ]; then
awk -v f_header=$f_header -v f_records=$f_records '{if ($0 ~ /^##/){print $0 > f_header} else {gsub("\\./\\.","0/0",$0);print $0 > f_records}}' ${f_addIDvcf}
mv ${f_addIDvcf} ${f_addIDvcf}_raw
cat $f_header $f_records > ${f_addIDvcf}
rm -f ${f_addIDvcf}_raw
else
	#only change ./. to 0/0 in f_addIDvcf
	found=$(grep "\./\." ${f_records}|head -1|wc -l)
	if [ $found -gt 0 ]; then
	mv $f_records ${f_records}_raw
	awk '{gsub("\\./\\.","0/0",$0);print $0}' ${f_records}_raw > ${f_records}
	fi
mv ${f_addIDvcf} ${f_addIDvcf}_raw
cat $f_header $f_records > ${f_addIDvcf}
rm -f ${f_addIDvcf}_raw ${f_records}_raw
fi
#---generate some useful intermediate files ed

#get allele frequency and reformat it
#	_ERVcaller_addID_freq.txt
pivcf=${cohort}_ERVcaller_addID.vcf
pofrq=${cohort}_ERVcaller_addID_freq
pi=${pofrq}.frq
po=${pofrq}.txt
vcftools --vcf $pivcf --freq --out $pofrq
awk 'NR==1{sub(/:/,"\t",$5);gsub(/[\{\}]/,"",$5);allel2=$5;sub(/ALLELE/,"alt_ALLELE",allel2);sub(/FREQ/,"alt_FREQ",allel2);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"allel2; exit}' $pi > $po
awk 'NR>1{sub(/:/,"\t",$5);gsub(/[\<\>]/,"",$6);split($6,arr,":");arrlen=length(arr);new6=arr[2];i=3;while (i<arrlen){new6=new6":"arr[i];i++};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"new6"\t"arr[arrlen]}' $pi >> $po


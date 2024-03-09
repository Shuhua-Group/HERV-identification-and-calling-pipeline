#!/usr/bin/bash
#cmd: bash sp01.sh PATH/TO/YOUR.s.config sampleID


configPath=$1
. $configPath
i=$2

cd $EC_out_dir

#make adpt.vcf
inDir=ori_vcf
IDprefix=""
IDsufffix=_adpt
awk 'BEGIN{i=1;j=1}{if(FNR==NR){{if(FNR>=2){split($15,a,/[^naDE]+/);arr[i]="_"a[1]"_"a[2];gsub("D",2,arr[i]);gsub("E",1,arr[i]);gsub("na",0,arr[i]);i++;next}}} else{if($1 !~ /^#/){split($8,a1,/;/);a1_2new=a1[2]arr[j];j++;a1[2]=a1_2new;new_col8="";for (m=1;m<=length(a1);m++){new_col8=new_col8";"a1[m]};sub(";","",new_col8);$8=new_col8;gsub(" ","\t",$0)};print $0}}' other_output/${i}_temp/${i}_ERV.output3 $inDir/${i}.vcf > ${i}_adpt.vcf

#make .bcf for merge
bcftools view $IDprefix$i$IDsufffix.vcf | bcftools sort -O b9 -o $IDprefix$i${IDsufffix}_sorted.bcf
bcftools index $IDprefix$i${IDsufffix}_sorted.bcf && \
mv ${i}_adpt\.vcf adpt_vcf/

cd $postEC_wd

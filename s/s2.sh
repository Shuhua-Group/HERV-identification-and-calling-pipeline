#!/usr/bin/bash


configPath=$1
. $configPath


#shell argument in
dataset=fltSetO1
infile_feature=${cohort}_ERVcaller_${dataset}

#out basenames
outfile_feature1=${infile_feature}_step1

outfile_feature3=${outfile_feature1}_represenID
outfile_feature5=${outfile_feature3}_clustermergedReGT

f_represenPosVcf=${outfile_feature5}.vcf




#make vcf namelist for merge
for i in reGT/*readSupportGTreplaced\.vcf; do 
	bcftools sort -O b9 -o ${i/vcf/bcf} $i
	bcftools index ${i/vcf/bcf}
done

f_reGTlist="reGT_bcf_list.txt"
ls reGT/*readSupportGTreplaced\.bcf > $f_reGTlist

#merge vcf with representative ID only.
bcftools merge --info-rules INFOR:join -l $f_reGTlist -O v -o ${f_represenPosVcf}
cp ${f_represenPosVcf} ${f_represenPosVcf}_temp
awk '{gsub("./.:.:.:.:.","./.:.:.,.,.:.:.",$0); print $0}' ${f_represenPosVcf}_temp > ${f_represenPosVcf}

vcftools --vcf ${f_represenPosVcf} --freq --out ${outfile_feature5}_frq
pi=${outfile_feature5}_frq.frq
po=${outfile_feature5}_frq.tsv
awk 'NR==1{sub(/:/,"\t",$5);gsub(/[\{\}]/,"",$5);allel2=$5;sub(/ALLELE/,"alt_ALLELE",allel2);sub(/FREQ/,"alt_FREQ",allel2);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"allel2; exit}' $pi > $po
awk 'NR>1{sub(/:/,"\t",$5);gsub(/[\<\>]/,"",$6);split($6,arr,":");arrlen=length(arr);new6=arr[1];i=2;while (i<arrlen){new6=new6":"arr[i];i++};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"new6"\t"arr[arrlen]}' $pi >> $po

vcftools --vcf ${f_represenPosVcf} --hardy --out ${outfile_feature5}_hwe
vcftools --missing-site --vcf ${f_represenPosVcf} --out ${outfile_feature5}_missing

bcftools view -H ${outfile_feature5}.vcf|awk '{print $3}' > ${outfile_feature5}_IDoriOrder.txt

rm $f_reGTlist || delete $f_reGTlist
rm reGT/*readSupportGTreplaced\.bcf* || delete reGT/*readSupportGTreplaced\.bcf*
rm ${f_represenPosVcf}_temp || delete ${f_represenPosVcf}_temp


##make GQ annotation after merge
#make_clustermergedReGT_GQannotation=${postERVcaller}/s/make_clustermergedReGT_GQannotation.py
#${python3} ${make_clustermergedReGT_GQannotation} ${outfile_feature5}_GQannotation.txt

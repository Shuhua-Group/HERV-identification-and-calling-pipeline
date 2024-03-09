import pandas_plink
import pandas as pd
import sys


in_plink=sys.argv[1]
out_pkl=sys.argv[2]

#reformat the plink output so that the allele 2 is the alternative allele but not the minor allele.
snp_info,sample_info,genotypes = pandas_plink.read_plink(in_plink)
genotype_mat = genotypes.compute()
df_geno=pd.DataFrame(data=genotype_mat,columns=sample_info['iid'],index=snp_info['snp'])

#genotypes are recorded as one in 2,1,0, where "2" stands for homozygous allele 1 (the reference allele), so change it to the widely-accepted formating 0,1,2.
df_geno=df_geno.applymap(lambda x: 2-x)

df_geno.to_pickle(out_pkl)
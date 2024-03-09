#filter vcf data set.
#will clear contents in INFO columns of the records in .vcf. This is simply to avoid potential format error warnings when the vcf file is used as input.


import pandas as pd
import sys, re


in_vcfRec=sys.argv[1]
in_pos_cluster=sys.argv[2]
in_clustermerged=sys.argv[3]
out_vcfRec=sys.argv[4]

in_posID_toDeleted=None
out_pos=None
out_posID=None
try:
    in_posID_toDeleted=sys.argv[5]
    out_pos=sys.argv[6]
    out_posID=sys.argv[7]
except:
    pass


if in_posID_toDeleted is not None:
    step=2
elif re.search("readSupportGT.vcf", in_vcfRec) is None:
    step=1
else:
    step=1.5


df_vcf=pd.read_csv(in_vcfRec,sep='\t',index_col=None)
df_vcf['ALT']=['<INS>']*len(df_vcf)


df_cluster=pd.read_csv(in_pos_cluster,sep='\t',index_col=None)
if len(df_vcf)>len(df_cluster):
    howmerge='right'
else:
    howmerge='left'
df_vcf=df_vcf.merge(df_cluster[['pos_id','cluster']],left_on='ID',right_on='pos_id',how=howmerge)
df_vcf=df_vcf.drop(columns=['pos_id'])

df_cluster_represen=pd.read_csv(in_clustermerged,sep='\t',index_col=None)
df_cluster_represen=df_cluster_represen.merge(df_vcf[['ID','REF']],left_on='pos_id',right_on='ID',how='left')
df_cluster_represen=df_cluster_represen.drop(columns=['ID'])
df_cluster_represen=df_cluster_represen.rename(columns={'pos':'represen_pos','pos_id':'represen_pos_id','REF':'represen_REF'})

if in_posID_toDeleted is not None:
    try:
        lis_posID_toDeleted=pd.read_csv(in_posID_toDeleted,sep='\t',header=None,index_col=None)[0].to_list()
    except:
        lis_posID_toDeleted=[]
    df_cluster_represen=df_cluster_represen.loc[~df_cluster_represen['represen_pos_id'].isin(lis_posID_toDeleted)]

if len(df_vcf)>len(df_cluster_represen):
    howmerge='right'
else:
    howmerge='left'
df_vcf=df_vcf.merge(df_cluster_represen[['cluster','represen_pos','represen_pos_id','represen_REF']],on='cluster',how=howmerge)


if step==1:
    df_vcf.iloc[:,:-4].to_csv(out_vcfRec0,header=1,index=None,sep='\t')
else:
    df_vcf_modi=df_vcf.copy()
    df_vcf_modi['POS']=df_vcf_modi['represen_pos'].astype(int)
    df_vcf_modi['ID']=df_vcf_modi['represen_pos_id']
    df_vcf_modi['REF']=df_vcf_modi['represen_REF']
    df_vcf_modi['INFO']=['.']*len(df_vcf_modi)
    df_vcf_modi.iloc[:,:-4].to_csv(out_vcfRec,header=1,index=None,sep='\t')
    df_vcf_modi[['#CHROM','POS','ID']].to_csv(out_pos,header=0,index=None,sep='\t')
    df_vcf_modi['ID'].to_csv(out_posID,header=0,index=None,sep='\t')


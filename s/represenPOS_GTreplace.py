#!/home/zhengwanjing/venv/bin/python3

#This is an independent script. Not called in other scripts.
#This is v3.1. added Line 39-40. changed definition of getBestRec(df_in,GT_1num), goodmap1vs0index
#This is v3.2. added flanking virtual sites to assist typing 1/0.

import pandas as pd
import os,sys
#import getopt

argv = sys.argv[1:]

if '-h' in argv:
    print('Version = v3. For per sample.\n\
    cmd format (4 arguments):\n\
    python3 represenPos_GTreplace.py \n\
    sample_ID \n\
    path_readSupportGT_vcf \n\
    path_step1_pos_cluster_info_txt \n\
    path_step1_clustermerged_cluster_info_txt\n\n\
    Will rewrite one output VCF to ./reGT/')
    quit()
#opts, args = getopt.getopt(argv, 'h:')
#if opts[0][0]=='-h':
#    print('cmd format:')
#    exit()




sampleID=sys.argv[1]

#intermediate file names p1
in_vcf_rec=sampleID+'_readSupportGT.vcf_rec'

#infile names p2
f_posCluster_info=sys.argv[3] #e.g. coh1_ERVcaller_filtSetO1_step1_pos_cluster_info.txt
f_clustermerged_info=sys.argv[4] #e.g. coh1_ERVcaller_filtSetO1_step1_clustermerged_cluster_info.txt

#out infile names
out_vcf_rec=in_vcf_rec.replace('readSupportGT','readSupportGTreplaced')
out_vcf=out_vcf_rec.replace('vcf_rec','vcf')
f_GQstrong_out=out_vcf_rec.replace('.vcf_rec','GQ_info.txt')

#make input per-sample vcf_record file (in_vcf_rec)
os.system("readSupportGTvcf="+sys.argv[2]+"; \
perSampleReadSupportGTvcf="+sampleID+"_readSupportGT.vcf; bcftools view -H -s "+sampleID+" -o $perSampleReadSupportGTvcf $readSupportGTvcf; \
mv $perSampleReadSupportGTvcf ${perSampleReadSupportGTvcf/.vcf/.vcf_rec}")


#reusable
lis_missingSymbols=['.','./.','1/.','0/.','./1','./0']#seems only './.' exists.

#import step1_pos-cluster-represenID info
df_posCluster_info=pd.read_csv(f_posCluster_info,index_col=None,sep='\t')
df_posCluster_info['chrpos']=df_posCluster_info[['chr','pos']].apply(
    lambda x: x[0]+':'+str(x[1]), axis=1)

df_clustermerged_info=pd.read_csv(f_clustermerged_info,index_col=None,sep='\t').iloc[:,:9]
df_clustermerged_info=df_clustermerged_info.rename(columns={'pos_id':'represenID'})

df_posCluster_info=df_posCluster_info.merge(
    df_clustermerged_info[['represenID','cluster','cluster_bg','cluster_ed']],
    on='cluster',how='left')#,'element'
lis_represenID=df_clustermerged_info['represenID'].tolist()

#make vcf df with pos-step1cluster-represenID info
vcfheader=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GT']
df_indiv_vcf=pd.read_csv(in_vcf_rec,names=vcfheader,sep='\t').sort_values(by=['CHROM','POS','ALT'])
df_indiv_vcf['pos_id']=df_indiv_vcf[['CHROM','POS','ALT']].apply(
    lambda x: x[0]+':'+str(x[1])+':'+x[2].split(':')[-1].replace('>',''), axis=1)
df_indiv_vcf['pos_id']=df_indiv_vcf[['pos_id','ID','ALT']].apply(lambda x: x[1] if x[2]=='INS' else x[0], axis=1)
#<--for virtualPos
df_indiv_vcf['chrpos']=df_indiv_vcf[['CHROM','POS']].apply(
    lambda x: x[0]+':'+str(x[1]), axis=1)

#Important note!! df_posCluster_info contain rows of same chrpos but different clusters (represenIDs) due to multiple elements/direction/... detected that produced a cluster boundary and were asigned to different cluster. df_clustermerged_info is after cluster merge and the represenPos are luckily no such cluster boundary ID were re-used for adjacent clusters, so doesn't have this issue.

df_indiv_vcf=df_indiv_vcf.merge(df_posCluster_info[['pos_id','cluster','represenID']],
                  on='pos_id', how='left')#,'element'
df_indiv_vcf=df_indiv_vcf.drop_duplicates()

#keep filtSet pos only
df_indiv_vcf=df_indiv_vcf.loc[~df_indiv_vcf['cluster'].isna()]
df_indiv_vcf['cluster']=df_indiv_vcf['cluster'].astype(int)

#add/change 4 info columns
#df_indiv_vcf['is_non_virtual']=df_indiv_vcf[['POS','cluster_bg','cluster_ed']].apply(lambda x: x[0]>=x[1] and x[0]<=x[2], axis=1)
df_indiv_vcf['pos_id']=df_indiv_vcf[['CHROM','POS','ALT','pos_id']]\
.apply(lambda x: x[0]+':'+str(x[1])+':'+x[2] if x[2]=='INS' else x[3], axis=1)
df_indiv_vcf['GT_short']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[0])
df_indiv_vcf['GT_short']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[0])
df_indiv_vcf['GQ']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[1]).apply(lambda x: -1 if x=='.' else int(x)).astype(int)


lis_repPosMissGT=df_indiv_vcf.loc[(df_indiv_vcf['GQ']==-1)&
                                 (df_indiv_vcf['pos_id'].isin(lis_represenID)),
                                 'represenID'].unique().tolist()


df_indiv_vcf_backup1=df_indiv_vcf.copy()
del df_posCluster_info,df_clustermerged_info




#make a initial dataframe for replacing missing GT of RepresenID
df_new_GTforRepresenID=pd.DataFrame()
df_indiv_vcf=df_indiv_vcf_backup1.copy()

df_indiv_vcf['GT_1num']=df_indiv_vcf['GT_short'].apply(lambda x: -1 if x in lis_missingSymbols else int(x.split('/')[0])+int(x.split('/')[1])).astype(int)

df_indiv_vcf['L00/L11']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[2].split(',') if x[2]!='.' else '.,.,.').apply(
    lambda x: round(float(x[0])/(float(x[2])+1e-6),5) if x[2]!='.' else -1)
df_indiv_vcf['L01/L11']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[2].split(',') if x[2]!='.' else '.,.,.').apply(
    lambda x: round(float(x[1])/(float(x[2])+1e-6),5) if x[2]!='.' else -1)

df_indiv_vcf['noIns_support']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[-2]).apply(lambda x: 0 if x=='.' else int(x)).astype(int)
df_indiv_vcf['ins_support']=df_indiv_vcf['GT'].apply(lambda x: x.split(':')[-1]).apply(lambda x: 0 if x=='.' else int(x)).astype(int)
df_indiv_vcf['readsSupport']=df_indiv_vcf[['noIns_support','ins_support']].sum(axis=1)

df_indiv_vcf['repPos']=df_indiv_vcf['represenID'].apply(lambda x: int(x.split(':')[1])).astype(int)
df_indiv_vcf['distance_to_repPos']=abs(df_indiv_vcf['POS']-df_indiv_vcf['repPos'])
df_indiv_vcf['distance_to_repPos_class']=df_indiv_vcf['distance_to_repPos']//20
df_indiv_vcf['DfamAcc_distance']=df_indiv_vcf[['pos_id','represenID']].apply(
    lambda x: abs(int(x[0].split('DF')[1].split('.')[0])-int(x[1].split('DF')[1].split('.')[0])) if len(x[0].split('DF'))>1 else None, axis=1)
df_indiv_vcf['GQstrong']=[0]*len(df_indiv_vcf)


#define functions useful for replacing missing GT of RepresenID
 #Get good records per GT
def getBestRec(df_in,GT_1num):
    df_work=df_in.loc[df_in['GT_1num']==GT_1num].copy()
    if len(df_work)==0:
        return GT_1num, None, None
    
    df_work1=df_work.loc[df_work['ALT']!='INS'].copy()#using the virtual sites
    df_work1=df_work1.sort_values(by=['GQ','distance_to_repPos_class','readsSupport','DfamAcc_distance'],
        ascending=[False,True,False,True],na_position='last')
    if len(df_work1)==0:
        return GT_1num, None, None
    
    df_candiRec0=df_work1.iloc[[0],:].copy()
    
    df_work=df_work.sort_values(by=['GQ','readsSupport','distance_to_repPos','DfamAcc_distance'],
        ascending=[False,False,False,True],na_position='last') #the more distant the better to show the environment near but not be the insertion
    df_topSupport=df_work.iloc[[0],:].copy()
    top_GQ=df_topSupport['GQ'].iat[0]
    if len(df_topSupport)>1:
        readsSupport_rubbed=int(round(df_work.loc[df_work['GQ']>top_GQ - 0.8 * (top_GQ/10) ** 2,'readsSupport'].median(),0))
        noIns_support_rubbed=int(round(df_topSupport['noIns_support'].iat[0]*readsSupport_rubbed/df_topSupport['readsSupport'].iat[0]))
        df_topSupport.loc[df_topSupport.index[0],'readsSupport']=readsSupport_rubbed
        df_topSupport.loc[df_topSupport.index[0],'noIns_support']=noIns_support_rubbed
        
    return GT_1num, df_candiRec0, df_topSupport


#do replacing missing GT of RepresenID
dic_GT1num_GQveryStrong={0:25,1:30,2:30}
dic_GT1num_GQstrong={0:20,1:20,2:20}
dic_GT1num_GQmoderate={0:15,1:20,2:20}
dic_GT1num_GQleftlim={0:10,1:10,2:10}

df_new_GTforRepresenID=pd.DataFrame()

for reppos in lis_represenID:
    df_temp0=df_indiv_vcf.loc[df_indiv_vcf['represenID']==reppos].copy()
    df_temp=df_temp0.loc[df_temp0['ALT']!='INS']#using the virtual sites
    if len(df_temp)==0:
        continue
    if len(df_temp)==1:
        df_candiRec=df_temp.copy()
        continue
    else:
        df_candiRec0=df_temp.loc[df_indiv_vcf['pos_id']==reppos].copy()
        df_candiRec=df_candiRec0.copy()
    
    #important: GT selection by sorting and booleans
    #important: GT selection by sorting and booleans
    #important: GT selection by sorting and booleans
    df_bestRecGTmiss=getBestRec(df_temp0,-1)[1]
    df_bestRecGT0,df_topSupportGT0=getBestRec(df_temp0,0)[1:]
    df_bestRecGT1=getBestRec(df_temp0,1)[1]
    df_bestRecGT2=getBestRec(df_temp0,2)[1]
    
    
    df_bestRecGT_non1=None
    changeCode=0
    if df_bestRecGT2 is not None:
        df_bestRecGT_non1=df_bestRecGT2
        df_candiRec=df_bestRecGT_non1.copy()
        changeCode+=1
        if df_bestRecGT0 is not None:
            goodmap2vs0index=df_bestRecGT2['readsSupport'].iat[0]/(df_topSupportGT0['readsSupport'].iat[0]+1e-6)
            if df_bestRecGT2['GQ'].iat[0] < dic_GT1num_GQleftlim[2] and goodmap2vs0index < 0.1:
                df_bestRecGT_non1=df_bestRecGT0.copy()
                df_candiRec=df_bestRecGT_non1.copy()
                changeCode+=10
    elif df_bestRecGT0 is not None:
        df_bestRecGT_non1=df_bestRecGT0.copy()
        df_candiRec=df_bestRecGT_non1.copy()
        changeCode+=100
    
    if df_bestRecGT1 is not None:
        good01map=False
        good01GQleftlimit=1e6
        good01suppl=True
        if df_bestRecGT_non1 is not None:
            if df_bestRecGT0 is not None and (df_bestRecGT_non1.iloc[0] == df_bestRecGT0.iloc[0]).all():
                goodmap1vs0index=(df_topSupportGT0['noIns_support'].iat[0]-df_bestRecGT1['noIns_support'].iat[0])/\
                                (df_topSupportGT0['noIns_support'].iat[0]+1e-6) *\
                                int(df_bestRecGT1['readsSupport'].iat[0]> 0.5 * df_topSupportGT0['readsSupport'].iat[0] and \
                                (df_bestRecGT1['GQ'].iat[0]>=dic_GT1num_GQveryStrong[1] or df_bestRecGT1['readsSupport'].iat[0]<=1.2 * df_topSupportGT0['readsSupport'].iat[0]))
                                #<--used evidence for: a value showing replaced(reduced) non-insertion supports to insertion(anti multi-mapping), informative enough number of reads, not too many reads(anti multi-mapping)
                if goodmap1vs0index>=0.2 or df_bestRecGT1['GQ'].iat[0]>=df_bestRecGT_non1['GQ'].iat[0]:
                    good01map=True
                good01GQleftlimit=dic_GT1num_GQmoderate[1]
                if goodmap1vs0index<=0 or df_bestRecGT1['GQ'].iat[0]<dic_GT1num_GQleftlim[1]:
                    good01suppl=False
            else:
                #fac1=max(4,df_bestRecGT_non1['noIns_support'].iat[0]//10)
                #fac2=min(4,df_bestRecGT_non1['GQ'].iat[0]//10)
                good01GQleftlimit=df_bestRecGT_non1['GQ'].iat[0]
                if df_bestRecGT_non1['noIns_support'].iat[0]>0 and df_bestRecGT_non1['L01/L11'].iat[0]>0.5:
                #if goodmap1vs2index>=0.2:
                    good01map=True                    
        if (df_bestRecGT1['GQ'].iat[0] > good01GQleftlimit or good01map) and good01suppl:
            df_candiRec=df_bestRecGT1.copy()
            changeCode+=1000
    
    df_new_GTforRepresenID=pd.concat([df_new_GTforRepresenID,df_candiRec],axis=0)


#do the replacement p1
df_new_GTforRepresenID['original_ID']=df_new_GTforRepresenID['pos_id']
df_new_GTforRepresenID['ID']=df_new_GTforRepresenID['represenID']
df_new_GTforRepresenID['pos_id']=df_new_GTforRepresenID['ID']
df_new_GTforRepresenID['original_pos']=df_new_GTforRepresenID['POS']
df_new_GTforRepresenID['POS']=df_new_GTforRepresenID['represenID'].apply(lambda x: int(x.split(':')[1]))
df_new_GTforRepresenID['POS']=df_new_GTforRepresenID['POS'].astype(int)

#make dictionary for represenPos reference and alt alleles in VCF for the replacement
df_represenID_REF_ALT=\
df_indiv_vcf.loc[(df_indiv_vcf['ALT']!='INS') & 
                 (df_indiv_vcf['ID'].isin(lis_represenID)),
                 ['ID','REF','ALT']].set_index('ID')
dic_represenID_REF_ALT=df_represenID_REF_ALT.to_dict(orient='index')
del df_represenID_REF_ALT

#do the replacement p2
df_new_GTforRepresenID['original_REF']=df_new_GTforRepresenID['REF']
df_new_GTforRepresenID['REF']=df_new_GTforRepresenID['represenID'].apply(lambda x: dic_represenID_REF_ALT[x]['REF'])
df_new_GTforRepresenID['original_ALT']=df_new_GTforRepresenID['ALT']
df_new_GTforRepresenID['ALT']=df_new_GTforRepresenID['represenID'].apply(lambda x: dic_represenID_REF_ALT[x]['ALT'])


#df_indiv_vcf_updated will contein the represenID positions and alelles only.
df_indiv_vcf_updated=df_indiv_vcf.loc[df_indiv_vcf['ALT']!='INS'].copy()
df_indiv_vcf_updated=\
df_indiv_vcf_updated.loc[(~df_indiv_vcf_updated['pos_id'].isin(df_new_GTforRepresenID['pos_id'])) | 
                         (df_indiv_vcf_updated.index.isin(df_new_GTforRepresenID.index))].copy()
#<--because the index for the same ID (value=pos_id) are different in the two input dfs, will make duplicates of ID in the updated df.
df_indiv_vcf_updated.update(df_new_GTforRepresenID)
df_indiv_vcf_updated=df_indiv_vcf_updated.loc[df_indiv_vcf_updated['ID'].isin(lis_represenID)]


df_indiv_vcf_updated['GQstrong']=\
df_indiv_vcf_updated[['GQ','GT_1num','GQstrong']].apply(
    lambda x: 1 if x[0]>=0 and x[0]>=dic_GT1num_GQstrong[x[1]] \
    else x[2],axis=1) 
df_indiv_vcf_updated['GQstrong']=\
df_indiv_vcf_updated[['GQ','GT_1num','GQstrong']].apply(
    lambda x: 2 if x[0]>=0 and x[0]>=dic_GT1num_GQveryStrong[x[1]] \
    else x[2],axis=1) 


for i in ('POS','cluster','GQ','GQstrong','GT_1num','noIns_support',
          'ins_support','readsSupport','repPos','distance_to_repPos',
          'distance_to_repPos_class'):
    df_indiv_vcf_updated[i]=df_indiv_vcf_updated[i].astype(int)





#export
os.system("\
if ! [ -d \"reGT\" ]; then\n\
	mkdir reGT\n\
fi\
")

df_vcf_out=df_indiv_vcf_updated.iloc[:,:10]
df_vcf_out.to_csv('reGT/'+out_vcf_rec,header=None,index=None,sep='\t')
df_GQstrong_out=df_indiv_vcf_updated[['ID','cluster','GQ','GQstrong','GT_1num']].rename(columns={'ID':'represenID'})
df_GQstrong_out.to_csv('reGT/'+f_GQstrong_out,header=1,index=None,sep='\t')

os.system("readSupportGTvcf="+sys.argv[2]+"; \
bcftools view -h -s "+sampleID+" -o reGT/"+sampleID+"_header.vcf $readSupportGTvcf; \
cat reGT/"+sampleID+"_header.vcf reGT/"+out_vcf_rec+" > reGT/"+out_vcf+"; \
rm reGT/"+sampleID+"_header.vcf || delete reGT/"+sampleID+"_header.vcf; \
rm "+in_vcf_rec+"|| delete "+in_vcf_rec)

#os.system("rm reGT/"+out_vcf_rec+"|| delete reGT/"+out_vcf_rec)


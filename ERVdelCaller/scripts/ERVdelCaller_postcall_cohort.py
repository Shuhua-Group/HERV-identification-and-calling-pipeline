import pandas as pd
import numpy as np
import sys,re,os


confp=sys.argv[1]
confd=re.sub(r'/*[^\s/]*\.config','',confp)
tpcp=confd+'/tpc.py'
os.system("cp "+confp+ " "+tpcp)
sys.path.append(confd)
from tpc import *


#os.system("bash "+ERVdelCaller+'/scripts/ERVdelCaller_merge.sh '+confp+' || echo -e \"\\nMerge failed!\\n\"')


data_dir_i=EdC_wd+'/merge_filter/'
data_dir_o=EdC_wd+'/'
coh=cohort


df_family=pd.read_csv(postERVcaller+'/p/element.info',sep='\t',header=0)
df_family=df_family.iloc[:,[0,1,2,3,4,5,6,8]]

lis_dfam_elememts=df_family['element'].tolist()

dic_element_accPartClassOrg=df_family.set_index('element')\
[['acc','ERV_part','312_class','organism']].apply(lambda x:list(x),axis=1)


def get_dfam_item(in_str):
    dfam_ele=True
    HERVdIntInFull='unknown'
    try:
        out=dic_element_accPartClassOrg[in_str]
        HERVdIntInFull='DfamPartOrFull'
        return in_str,out,HERVdIntInFull
    except:
        dfam_ele=False
    
    for i in ('-int','_int','_'):
        if len(i)>1:
            in_str_new=in_str.replace(i,'')
            if in_str_new != in_str:
                HERVdIntInFull='INT'
                dic_HERVdFull_HERVdInt[in_str_new]=in_str
        else:
            in_str_new=in_str.rstrip(i)
        if in_str_new in lis_dfam_elememts:
            out=dic_element_accPartClassOrg[in_str_new]
            dfam_ele=True
            break
        else:
            in_str_new=in_str+i
            if in_str_new in lis_dfam_elememts:
                HERVdIntInFull='FULL'
                dic_HERVdFull_HERVdInt[in_str]=in_str_new
                acc=dic_element_accPartClassOrg[in_str_new]
                dfam_ele=True
                break

    if dfam_ele==True:
        return in_str_new,out,HERVdIntInFull
    else:
        return '',['','','',''],HERVdIntInFull


#ERV deletion table
df_del_info0=pd.read_csv(data_dir_i+coh+'_ERVdelSV_info.txt',sep='\t')
dic_HERVdFull_HERVdInt={}

df_del_info0['ori_index']=df_del_info0.index.values
df_del_info0['dist_abs_sum']=df_del_info0[['left_distance','right_distance']]\
.apply(lambda x: abs(x[0])+abs(x[1]),axis=1)

df_del_info=df_del_info0.copy()

#If a new column name does not include "HERVd", it refers to Dfam.
df_del_info['HERVd_elements_list']=\
df_del_info['reference_interval_ID'].apply(lambda x: re.sub(r'^ERV_\d+_','',x).split(','))

df_del_info['elements']=df_del_info['HERVd_elements_list']\
.apply(lambda x: sorted(list(set([get_dfam_item(x[i])[0] for i in range(len(x))]))))
df_del_info['312_classes']=df_del_info['HERVd_elements_list']\
.apply(lambda x: sorted(list(set([get_dfam_item(x[i])[1][2] for i in range(len(x))]))))
df_del_info['ERV_parts_list']=df_del_info['HERVd_elements_list']\
.apply(lambda x: [get_dfam_item(x[i])[1][1] for i in range(len(x))])
df_del_info['ERV_part']=df_del_info['HERVd_elements_list']\
.apply(lambda x: sorted(list(set([get_dfam_item(x[i])[1][1] for i in range(len(x))]))))

df_del_info['element']=df_del_info['HERVd_elements_list']\
.apply(lambda x: ','.join(sorted(list(set([get_dfam_item(x[i])[0] for i in range(len(x))])))))
df_del_info['312_class']=df_del_info['HERVd_elements_list']\
.apply(lambda x: ','.join(sorted(list(set([get_dfam_item(x[i])[1][2] for i in range(len(x))])))))
df_del_info['ERV_part']=df_del_info['HERVd_elements_list']\
.apply(lambda x: '-'.join(sorted(list(set([get_dfam_item(x[i])[1][1] for i in range(len(x))])),reverse=True)))
df_temp=df_del_info[['ERV_parts_list','HERVd_elements_list']].copy()
df_temp['HERVd_parts_list']=df_temp['HERVd_elements_list']\
.apply(lambda x: [get_dfam_item(x[i])[2] for i in range(len(x))])
df_temp['HERVd_parts_list_v2']=df_temp.apply(
    lambda x: ['FULL' if x['HERVd_elements_list'][i] in dic_HERVdFull_HERVdInt.keys() 
               else x['HERVd_parts_list'][i] for i in range(len(x['HERVd_elements_list']))],
    axis=1)
df_temp['HERVd_parts_list_v3']=df_temp.apply(
    lambda x: [x['ERV_parts_list'][i] if x['HERVd_parts_list_v2'][i] == 'DfamPartOrFull' 
               else x['HERVd_parts_list_v2'][i] for i in range(len(x['HERVd_parts_list_v2']))],
    axis=1)
df_del_info['HERVd_part']=df_temp['HERVd_parts_list_v3']\
.apply(lambda x: '-'.join(sorted(list(set(x)),reverse=True)))

lastrowidx=df_del_info.index[0]
max_dist_current_SV=df_del_info['dist_abs_sum'].iat[0]
df_del_info['dist_abs_sum_ismax']=[False]*len(df_del_info)
for rec in df_del_info.iloc[1:].iterrows():
    rowidx=rec[0]
    current_dist=df_del_info.loc[rowidx,'dist_abs_sum']
    if df_del_info.loc[rowidx,'SV_ID'] == df_del_info.loc[lastrowidx,'SV_ID']:
        if current_dist > max_dist_current_SV:
            max_dist_current_SV=current_dist
    else:
        max_dist_current_SV=current_dist
        df_del_info.loc[rowidx,'dist_abs_sum_ismax']=True
    lastrowidx=rowidx
df_del_info=df_del_info.loc[df_del_info['dist_abs_sum_ismax']]


#filter
df_SVcount_per_refERV=df_del_info[['reference_interval_ID','SV_ID']].\
groupby('reference_interval_ID').count().rename(columns={'SV_ID':'SV count'})
lis_multi_allelic_ERVs=df_SVcount_per_refERV.loc[df_SVcount_per_refERV['SV count']>1].index.tolist()
df_del_finflt_info=df_del_info.loc[(~df_del_info['reference_interval_ID'].isin(lis_multi_allelic_ERVs))]

#export
path_in_vcf=data_dir_i+coh+'_ERVdelSV.vcf'
path_out_finflt_id=data_dir_o+coh+'_ERVdel_finflt_SVid.txt'
path_out_finflt_vcf=data_dir_o+coh+'_ERVdelSV_finflt.vcf'

df_out=df_del_finflt_info[['chromosome', 'SV_ID', 'SV_start_position', 'SV_length',
       'reference_interval_ID', 'reference_interval_start_position',
       'reference_interval_length', 'left_distance', 'right_distance',
       'alt_freq','312_class','element','HERVd_part']]
df_out.to_csv(data_dir_o+coh+'_ERVdel_finflt_info.txt',header=1,index=0,sep='\t')
df_del_finflt_info['SV_ID'].to_csv(path_out_finflt_id,header=0,index=0,sep='\t')
os.system('bcftools filter -i ID=@'+path_out_finflt_id+' -o '+path_out_finflt_vcf+' '+path_in_vcf)


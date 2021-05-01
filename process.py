import pandas as pd
import numpy as np
import os
import re

#保存蛋白质的分组信息
protein_cluster = {}
with open('protein.clstr','r')as f:
    for line in f:
        if line.startswith('>'):
            cluster = line[1:].rstrip()
        else:
            index = line.split('>')[1].find('.')
            protein = line.split('>')[1][:index]
            protein_cluster[protein] = cluster

#保存RNA的分组信息
RNA_cluster = {}
cluster_RNA = {}
with open('RNA.clstr','r')as f:
    for line in f:
        if line.startswith('>'):
            cluster = line[1:].rstrip()
        elif line[0].isdigit() :
            RNA,dataset = line.split('>')[1].split('_')
            dataset = dataset[:dataset.find('.')]
            RNA_cluster[RNA+"_"+dataset] = cluster
            #print(RNA+" "+cluster+" "+dataset)

#获取每个数据集出现在每个cluster_pair中的次数
def get_occuranceTimes(dataset_frame,dataset_name,final_interactions):
    for index,row in dataset_frame.iterrows():
        if(row['RNA'] not in ['YBL039W-A', 'YBL101W-A', 'YFL057C', 'YIR044C', 'YAR062W', 'YNL097C-A']):
            cluster_pair =(RNA_cluster[str(row['RNA'])+"_"+dataset_name],protein_cluster[str(row['protein'])+"_"+dataset_name])
        if final_interactions.__contains__(cluster_pair):
            final_interactions[cluster_pair][dataset_name]+=1
        else:
            final_interactions[cluster_pair] = {'RPI488':0,'RPI369':0,'RPI2241':0,'RPI7317':0,'RPI13254':0,
                                                'RPI1807':0,'RPI1446':0,'NPInter10412':0,'ZeaMays':0,'ArabidopsisThaliana':0
                                               }
            final_interactions[cluster_pair][dataset_name]+=1
    return final_interactions

final_interactions = {}
RPI488 = pd.read_table('RPI488_all.txt')
RPI488 = RPI488[RPI488['label']==1]
final_interactions1 = get_occuranceTimes(RPI488,'RPI488',final_interactions)

RPI369 = pd.read_table('RPI369_all.txt',header= None,names=['protein','RNA','label'])
RPI369 = RPI369[RPI369['label']==1]
final_interactions2 = get_occuranceTimes(RPI369,'RPI369',final_interactions1)

RPI2241 = pd.read_table('RPI2241_pairs.txt',header=None,names = ['protein','RNA','label'])
RPI2241 = RPI2241[RPI2241['label']==1]
final_interactions3 = get_occuranceTimes(RPI2241,'RPI2241',final_interactions2)

RPI1807 = pd.read_table('RPI1807_PositivePairs.csv')
final_interactions4 = get_occuranceTimes(RPI1807,'RPI1807',final_interactions3)

RPI1446 = pd.read_table('RPI1446_pairs.txt',header=None,names=['protein','RNA','label'])
RPI1446 = RPI1446[RPI1446['label']==1]
final_interactions5 = get_occuranceTimes(RPI1446,'RPI1446',final_interactions4)

RPI13254 = pd.read_table('RPI13254_pairs.txt')
final_interactions6 = get_occuranceTimes(RPI13254,'RPI13254',final_interactions5)

NPInter10412 = pd.read_table('NPInter10412_dataset.txt')
final_interactions7 = get_occuranceTimes(NPInter10412,'NPInter10412',final_interactions6)

RPI7317 = pd.read_csv('RPI7317.csv')
final_interactions8 = get_occuranceTimes(RPI7317,'RPI7317',final_interactions7)

ZeaMays = pd.read_table('Zea mays_all.txt')
final_interactions9 = get_occuranceTimes(ZeaMays,'ZeaMays',final_interactions8)

ArabidopsisThaliana = pd.read_table('Arabidopsis thaliana_all.txt')
final_interactions10 = get_occuranceTimes(ArabidopsisThaliana,'ArabidopsisThaliana',final_interactions9)
#print(final_interactions10)


with open('interaction_table.csv','a')as sheet:
    sheet.write(",RPI488,RPI369,RPI2241,RPI7317,RPI13254,RPI1807,RPI1446,NPInter10412,ZeaMays,ArabidopsisThaliana\n")
    for dataset1 in ['RPI488','RPI369','RPI2241','RPI7317','RPI13254',
                    'RPI1807','RPI1446','NPInter10412','ZeaMays','ArabidopsisThaliana']:
        sheet.write(dataset1+',')
        for dataset2 in ['RPI488','RPI369','RPI2241','RPI7317','RPI13254',
                    'RPI1807','RPI1446','NPInter10412','ZeaMays','ArabidopsisThaliana']:
            cnt = 0
            for key,value in final_interactions10.items():
                cnt+=min(value[dataset1],value[dataset2])
            sheet.write(str(cnt)+',')
        sheet.write('\n')

    sheet.close()


print(len(set(final_interactions10.keys())))

import sys,re,os

contron = sys.argv[1]
treatment = sys.argv[2]
rep = sys.argv[3]
regulate_dir = sys.argv[4]

gene = {}
title = ['GeneName','BriefDescription','GoAcession','TranscriptionFactor','ath_orth_gene']

contron_number = contron.split('/')
treatment_number = treatment.split('/')

TranscriptionFactor_gene = {}
with open(sys.argv[5],'r') as TranscriptionFactor:
    for i in TranscriptionFactor:
        ent = i.strip().split('\t')
        if ent[0].split('.')[0] not in TranscriptionFactor_gene:
            TranscriptionFactor_gene[ent[0].split('.')[0]] = ent[1]

ath_ort_gene = {}
with open(sys.argv[6],'r') as ath_orh:
    for i in ath_orh:
        ent = i.strip().split('\t')
        ath_ort_gene[ent[0].split('.')[0]] = ent[1]

for ii in range(len(contron_number)):
    for jj in range(len(treatment_number)):
        for kk in range(1,int(rep)+1):
            with open(regulate_dir+'/'+contron_number[ii]+'_vs_'+treatment_number[jj]+'_'+str(kk)+'.txt') as regulate:
                for i in regulate:
                    if i.startswith('GeneName'):
                        continue
                    else:
                        ent = i.strip().split('\t')
                        if ent[0] not in gene:
                            if ent[0] in TranscriptionFactor_gene:
                                if ent[0] in ath_ort_gene:
                                    gene[ent[0]] = ent[-3]+'\t'+ent[-2]+'\t'+TranscriptionFactor_gene[ent[0]]+'\t'+ath_ort_gene[ent[0]]
                                else:
                                    gene[ent[0]] = ent[-3]+'\t'+ent[-2]+'\t'+'nodcoding'+'\t'+'nodcoding'
                            else:
                                if ent[0] in ath_ort_gene:
                                    gene[ent[0]] = ent[-3]+'\t'+ent[-2]+'\t'+'NotTF'+'\t'+ath_ort_gene[ent[0]]
                                else:
                                    gene[ent[0]] = ent[-3]+'\t'+ent[-2]+'\t'+'NotTF'+'\t'+'nodcoding'


for ii in range(len(contron_number)):
    for jj in range(len(treatment_number)):
        for kk in range(1,int(rep)+1):
            with open(regulate_dir+'/'+contron_number[ii]+'_vs_'+treatment_number[jj]+'_'+str(kk)+'.txt') as regulate:
                title.append(contron_number[ii]+'_vs_'+treatment_number[jj]+'_'+str(kk))
                gene_in_this_rep = []
                for i in regulate:
                    if i.startswith('GeneName'):
                        continue
                    else:
                        ent = i.strip().split('\t')
                        if ent[0] in gene:
                            gene_in_this_rep.append(ent[0])
                            if float(ent[1]) < 1 and float(ent[2]) < 1 and float(ent[3]) < 1 and float(ent[4]) < 1 and float(ent[5]) <1 and float(ent[6]) < 1:
                                gene[ent[0]] += '\t'+ent[-1]+'/'+ent[-6]+'/FPKM<1'
                            else:
                                gene[ent[0]] += '\t'+ent[-1]+'/'+ent[-6]
                other_gene = list(set(gene.keys()).difference(set(gene_in_this_rep)))
                for j in other_gene:
                    gene[j] += '\t'+'------'

print ('\t'.join(title))
for i in gene:
    print (i+'\t'+gene[i])

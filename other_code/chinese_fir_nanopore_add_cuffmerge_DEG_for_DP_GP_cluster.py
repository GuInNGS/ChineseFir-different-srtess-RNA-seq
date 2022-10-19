import sys,re,os

DEG_gene = {}
with open(sys.argv[1],'r') as DEG:
    for i in DEG:
        if i.startswith('GeneName'):
            continue
        ent = i.strip().split('\t')
        if ent[0] not in DEG_gene:
            DEG_gene[ent[0]] = ''

with open(sys.argv[2],'r') as all_gene:
    for i in all_gene:
        if i.startswith('locus'):
            print (i.strip())
            continue
        ent = i.strip().split('\t')
        if ent[0] in DEG_gene:
            print (i.strip())


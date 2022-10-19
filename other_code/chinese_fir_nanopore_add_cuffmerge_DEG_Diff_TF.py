import sys,re,os

DEG_gene = {}
gene_family = sys.argv[3]

with open(sys.argv[1],'r') as DEG:
    for i in DEG:
        ent = i.strip().split('\t')
        if i.startswith('GeneName'):
            continue
        if re.search(gene_family,ent[3]):
            DEG_gene[ent[0]] = ''

with open(sys.argv[2],'r') as FPKM:
    for i in FPKM:
        if i.strip().startswith('locus'):
            print (i.strip())
        ent = i.strip().split('\t')
        if ent[0] in DEG_gene:
            print (i.strip())

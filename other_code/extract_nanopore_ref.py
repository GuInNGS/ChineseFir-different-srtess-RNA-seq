import sys,re,os

with open(sys.argv[1],'r') as gff:
    gff3 = gff.readlines()

with open(sys.argv[2],'r') as seq:
    ref_seq = seq.readlines()

gene_seq = {}
for i in range(0,len(ref_seq),2):
    ent = ref_seq[i].strip().split()[0].replace('>','')
    gene_seq[ent] = ref_seq[i+1].upper()

for i in range(len(gff3)):
    ent = gff3[i].strip().split()
    if len(ent) <= 2:
        continue
    if ent[2] == 'gene':
        gene_name = ent[-1].split('=')[-1]
        if gene_name in gene_seq:
            print ('>'+gene_name+'\n'+gene_seq[gene_name].strip())

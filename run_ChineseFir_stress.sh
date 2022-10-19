#source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate nanopore
#export PYTHONPATH=~/bigdata/hitseq/why/bin/tool/anaconda3/envs/nanopore/bin/
Genome_bowtie2_index=/home/lfgu/bigdata/Genome/chinaFir/Genome/bowtie2/Cunninghamia_lanceolata
Genome=/home/lfgu/bigdata/Genome/chinaFir/Genome/Cunninghamia_lanceolata_V1.fasta
Genome_fai=/home/lfgu/bigdata/Genome/chinaFir/Genome/Cunninghamia_lanceolata_V1.fasta.fai
gmap=/home/lfgu/bigdata/Genome/chinaFir/Genome/gmap
GFF=/home/lfgu/bigdata/Genome/chinaFir/gff/Cunninghamia_lanceolata_V1.gff3

input=/home/lfgu/bigdata/hitseq/SDX_fir/input
output=/home/lfgu/bigdata/hitseq/SDX_fir/output_nanopore
bin=/home/lfgu/bigdata/hitseq/why/bin/

Control_R1_rep1=$input/RNA-seq-rep1/C-SDX_L1_369369.R1.clean.fastq
Control_R2_rep1=$input/RNA-seq-rep1/C-SDX_L1_369369.R2.clean.fastq
Control_R1_rep2=$input/RNA-seq-rep2/C-SDX_L2_369369.R1.clean.fastq
Control_R2_rep2=$input/RNA-seq-rep2/C-SDX_L2_369369.R2.clean.fastq
Control_R1_rep3=$input/RNA-seq-rep3/C-SDX_L3_369369.R1.clean.fastq
Control_R2_rep3=$input/RNA-seq-rep3/C-SDX_L3_369369.R2.clean.fastq

Cold_1_R1_rep1=$input/RNA-seq-rep1/Cold-SDX-1_L1_370370.R1.clean.fastq
Cold_1_R2_rep1=$input/RNA-seq-rep1/Cold-SDX-1_L1_370370.R2.clean.fastq
Cold_1_R1_rep2=$input/RNA-seq-rep2/Cold-SDX-1_L2_370370.R1.clean.fastq
Cold_1_R2_rep2=$input/RNA-seq-rep2/Cold-SDX-1_L2_370370.R2.clean.fastq
Cold_1_R1_rep3=$input/RNA-seq-rep3/Cold-SDX-1_L3_370370.R1.clean.fastq
Cold_1_R2_rep3=$input/RNA-seq-rep3/Cold-SDX-1_L3_370370.R2.clean.fastq

Cold_2_R1_rep1=$input/RNA-seq-rep1/Cold-SDX-2_L1_371371.R1.clean.fastq
Cold_2_R2_rep1=$input/RNA-seq-rep1/Cold-SDX-2_L1_371371.R2.clean.fastq
Cold_2_R1_rep2=$input/RNA-seq-rep2/Cold-SDX-2_L2_371371.R1.clean.fastq
Cold_2_R2_rep2=$input/RNA-seq-rep2/Cold-SDX-2_L2_371371.R2.clean.fastq
Cold_2_R1_rep3=$input/RNA-seq-rep3/Cold-SDX-2_L3_371371.R1.clean.fastq
Cold_2_R2_rep3=$input/RNA-seq-rep3/Cold-SDX-2_L3_371371.R2.clean.fastq

Cold_3_R1_rep1=$input/RNA-seq-rep1/Cold-SDX-3_L1_372372.R1.clean.fastq
Cold_3_R2_rep1=$input/RNA-seq-rep1/Cold-SDX-3_L1_372372.R2.clean.fastq
Cold_3_R1_rep2=$input/RNA-seq-rep2/Cold-SDX-3_L2_372372.R1.clean.fastq
Cold_3_R2_rep2=$input/RNA-seq-rep2/Cold-SDX-3_L2_372372.R2.clean.fastq
Cold_3_R1_rep3=$input/RNA-seq-rep3/Cold-SDX-3_L3_372372.R1.clean.fastq
Cold_3_R2_rep3=$input/RNA-seq-rep3/Cold-SDX-3_L3_372372.R2.clean.fastq

PEG_1_R1_rep1=$input/RNA-seq-rep1/PEG-SDX-1_L1_374374.R1.clean.fastq
PEG_1_R2_rep1=$input/RNA-seq-rep1/PEG-SDX-1_L1_374374.R2.clean.fastq
PEG_1_R1_rep2=$input/RNA-seq-rep2/PEG-SDX-1_L2_374374.R1.clean.fastq
PEG_1_R2_rep2=$input/RNA-seq-rep2/PEG-SDX-1_L2_374374.R2.clean.fastq
PEG_1_R1_rep3=$input/RNA-seq-rep3/PEG-SDX-1_L3_374374.R1.clean.fastq
PEG_1_R2_rep3=$input/RNA-seq-rep3/PEG-SDX-1_L3_374374.R2.clean.fastq

PEG_2_R1_rep1=$input/RNA-seq-rep1/PEG-SDX-2_L1_375375.R1.clean.fastq
PEG_2_R2_rep1=$input/RNA-seq-rep1/PEG-SDX-2_L1_375375.R2.clean.fastq
PEG_2_R1_rep2=$input/RNA-seq-rep2/PEG-SDX-2_L2_375375.R1.clean.fastq
PEG_2_R2_rep2=$input/RNA-seq-rep2/PEG-SDX-2_L2_375375.R2.clean.fastq
PEG_2_R1_rep3=$input/RNA-seq-rep3/PEG-SDX-2_L3_375375.R1.clean.fastq
PEG_2_R2_rep3=$input/RNA-seq-rep3/PEG-SDX-2_L3_375375.R2.clean.fastq

PEG_3_R1_rep1=$input/RNA-seq-rep1/PEG-SDX-3_L1_376376.R1.clean.fastq
PEG_3_R2_rep1=$input/RNA-seq-rep1/PEG-SDX-3_L1_376376.R2.clean.fastq
PEG_3_R1_rep2=$input/RNA-seq-rep2/PEG-SDX-3_L2_376376.R1.clean.fastq
PEG_3_R2_rep2=$input/RNA-seq-rep2/PEG-SDX-3_L2_376376.R2.clean.fastq
PEG_3_R1_rep3=$input/RNA-seq-rep3/PEG-SDX-3_L3_376376.R1.clean.fastq
PEG_3_R2_rep3=$input/RNA-seq-rep3/PEG-SDX-3_L3_376376.R2.clean.fastq

SA_1_R1_rep1=$input/RNA-seq-rep1/SA-SDX-1_L1_377377.R1.clean.fastq
SA_1_R2_rep1=$input/RNA-seq-rep1/SA-SDX-1_L1_377377.R2.clean.fastq
SA_1_R1_rep2=$input/RNA-seq-rep2/SA-SDX-1_L2_377377.R1.clean.fastq
SA_1_R2_rep2=$input/RNA-seq-rep2/SA-SDX-1_L2_377377.R2.clean.fastq
SA_1_R1_rep3=$input/RNA-seq-rep3/SA-SDX-1_L3_377377.R1.clean.fastq
SA_1_R2_rep3=$input/RNA-seq-rep3/SA-SDX-1_L3_377377.R2.clean.fastq

SA_2_R1_rep1=$input/RNA-seq-rep1/SA-SDX-2_L1_378378.R1.clean.fastq
SA_2_R2_rep1=$input/RNA-seq-rep1/SA-SDX-2_L1_378378.R2.clean.fastq
SA_2_R1_rep2=$input/RNA-seq-rep2/SA-SDX-2_L2_378378.R1.clean.fastq
SA_2_R2_rep2=$input/RNA-seq-rep2/SA-SDX-2_L2_378378.R2.clean.fastq
SA_2_R1_rep3=$input/RNA-seq-rep3/SA-SDX-2_L3_378378.R1.clean.fastq
SA_2_R2_rep3=$input/RNA-seq-rep3/SA-SDX-2_L3_378378.R2.clean.fastq

SA_3_R1_rep1=$input/RNA-seq-rep1/SA-SDX-3_L1_379379.R1.clean.fastq
SA_3_R2_rep1=$input/RNA-seq-rep1/SA-SDX-3_L1_379379.R2.clean.fastq
SA_3_R1_rep2=$input/RNA-seq-rep2/SA-SDX-3_L2_379379.R1.clean.fastq
SA_3_R2_rep2=$input/RNA-seq-rep2/SA-SDX-3_L2_379379.R2.clean.fastq
SA_3_R1_rep3=$input/RNA-seq-rep3/SA-SDX-3_L3_379379.R1.clean.fastq
SA_3_R2_rep3=$input/RNA-seq-rep3/SA-SDX-3_L3_379379.R2.clean.fastq

if [ ! -f $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fa ]; then
   echo '---------------------------'
   echo "[`date`]  fastq2fasta"
   cat $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/fastq_pass/* > $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fq
   python $bin/python/nanopore_fastq_fasta_U2T.py $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fq $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fa
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----fastq2fasta-----'
fi

if [ ! -d "$output/caun" ]; then
   echo '---------------------------'
   echo "[`date`]  start"
   mkdir -p $output/caun
   canu useGrid=false -correct gnuplotImageFormat=svg corOutCoverage=10000 mhapSensitivity=high corMinCoverage=0 correctedErrorRate=0.16 overlapper=minimap minReadLength=200 stopOnLowCoverage=0.01  corThreads=30 cormmapThreads=30 corovlThreads=30  minOverlapLength=100 genomeSize=10g -p cla_nanopore -d $output/caun -nanopore-raw $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fa batMemory=250 batThreads=40
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----caun correct-----'
fi

if [ ! -d "$output/porechop" ]; then
   echo '---------------------------'
   echo "[`date`]  start"
   mkdir -p $output/porechop
   porechop --format fasta -t 47 -i $output/caun/cla_nanopore.correctedReads.fasta -o $output/porechop/cla.choped.fasta
   echo '---------------------------'
   echo "[`date`]  finsh"
   echo '-----------porechop--------'
fi

if [ ! -d "$output/lordec-correct_after_caun" ]; then
   echo '--------------------'
   echo '-----lordec-correct after caun-----'
   echo "[`date`]  start"
   mkdir -p $output/lordec-correct_after_caun
   #cat /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Cold_1_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Cold_2_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Cold_3_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Control_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/PEG_1_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/PEG_2_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/PEG_3_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/SA_1_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/SA_2_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/SA_3_rep1.fastq >/home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/lordec-correct.fastq
   #T=threads k=kmer_len(big genome=21,small genome=19) s=solid_threshold(2 or 3)  S=out statistics file  i=long read FASTA/Q file  2=short read FASTA/Q file(s)  o=corrected_read_file
   lordec-correct -T 35 -k 21 -s 3 -S $output/lordec-correct_after_caun/statistics.log -i $output/porechop/cla.choped.fasta -2 /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/lordec-correct.fastq -o $output/lordec-correct_after_caun/correct.fasta &
   wait
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----lordec-correct after caun-----'
fi

if [ ! -d "$output/minimap2_lordec-correct_after_caun" ]; then
   echo '--------------------'
   echo '-----minimap after lordec-correct after caun-----'
   echo "[`date`]  start"
   mkdir -p $output/minimap2_lordec-correct_after_caun
   minimap2 -ax splice -uf -k14 -t 35 $Genome $output/lordec-correct_after_caun/correct.fasta > $output/minimap2_lordec-correct_after_caun/aln.sam
   conda deactivate
   #
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate python27
   export PYTHONPATH=/home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/envs/python27/bin/python
   samtools view -bT $Genome -o $output/minimap2_lordec-correct_after_caun/aln.bam -@ 30 $output/minimap2_lordec-correct_after_caun/aln.sam
   samtools sort -@ 30 $output/minimap2_lordec-correct_after_caun/aln.bam >$output/minimap2_lordec-correct_after_caun/aln.sort.bam
   samtools index $output/minimap2_lordec-correct_after_caun/aln.sort
   conda deactivate
   #bam to gff
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate nanopore
   /home/lfgu/bigdata/hitseq/why/bin/tool/pinfish/spliced_bam2gff/spliced_bam2gff -t 40 $output/minimap2_lordec-correct_after_caun/aln.sort.bam >$output/minimap2_lordec-correct_after_caun/aln.gff
   conda deactivate
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----minimap2 lordec-correct after caun-----'
fi


if [ ! -d "$output/gmap_after_caun" ]; then
   echo '-----------------------'
   echo '---------gmap after_caun----------'
   echo "[`date`]  start"
   mkdir -p $output/gmap_after_caun
   #/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 25 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/lordec-correct_after_caun_cd_hit_est.fasta >$output/gmap_after_caun/gmap_after_caun.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xaa >$output/gmap_after_caun/xaa.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xab >$output/gmap_after_caun/xab.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xac >$output/gmap_after_caun/xac.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xad >$output/gmap_after_caun/xad.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xae >$output/gmap_after_caun/xae.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xaf >$output/gmap_after_caun/xaf.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xag >$output/gmap_after_caun/xag.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xah >$output/gmap_after_caun/xah.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xai >$output/gmap_after_caun/xai.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xaj >$output/gmap_after_caun/xaj.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xak >$output/gmap_after_caun/xak.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xal >$output/gmap_after_caun/xal.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xam >$output/gmap_after_caun/xam.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xan >$output/gmap_after_caun/xan.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xao >$output/gmap_after_caun/xao.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xap >$output/gmap_after_caun/xap.gff3
/home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/CD-HIT-EST/xaq >$output/gmap_after_caun/xaq.gff3
cat *.gff3 >$output/gmap_after_caun/all.gff3
rm $output/gmap_after_caun/x*
##fiter short read
python /home/lfgu/bigdata/hitseq/why/bin/python/fir_nanopore_longest.py  $output/gmap_after_caun/all.gff3 > $output/gmap_after_caun/uniq_gmap_after_caun.gff3
##merge SGS RNA-seq cuffmerge with nanopore
mkdir -p $output/cuffmerge_add_nanopore
cat /home/lfgu/bigdata/hitseq/SDX_fir/output/cuffmerge/merged.gff3 $output/gmap_after_caun/all.gff3 > $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore.gff3
python /home/lfgu/bigdata/hitseq/why/bin/python/fir_nanopore_longest.py $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore.gff3 >$output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.gff3
###extract seq for SGS RNA-seq cuffmerge with nanopore(nanopore seq are from Genome)
python /home/lfgu/bigdata/hitseq/why/bin/python/fir_SDX_extract_exon_gff3.py $Genome  $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.gff3 >$output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas
bowtie2-build $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas  $output/cuffmerge_add_nanopore/RSEM/index/nanopore
###extract seq for SGS RNA-seq cuffmerge with nanopore(nanopore seqs are real nanopore seq)
python /home/lfgu/bigdata/hitseq/why/bin/python/fir_SDX_extract_cufflink_exon_nanopore.py $Genome  $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.gff3 $output/nanopore_ref_seq/nanopore_ref.fa >$output/cuffmerge_add_nanopore/cuffmerge_add_real_nanopore_uniq.fas
mkdir -p $output/cuffmerge_add_nanopore/TransDecoder
/home/lfgu/tool/Trinity/trinityrnaseq_r20140413p1/trinity-plugins/transdecoder/TransDecoder -t $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas  -S --reuse --workdir $output/cuffmerge_add_nanopore
mv cuffmerge_add_nanopore_uniq.fas.transdecoder.* $output/cuffmerge_add_nanopore/TransDecoder
#nanopore overlap with gff
   python  /home/lfgu/bigdata/hitseq/why/bin/python/nanopore_gff_overlap.py $GFF $output/gmap_after_caun/culster_uniq_jbrowse.gff > $output/gmap_after_caun/culster_uniq_overlap_gff   
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----gmap_after_caun-----'
fi

if [ ! -d "$output/nanopore_ref_seq" ]; then
   echo '-----------------------'
   echo '---------nanopore seq----------'
   echo "[`date`]  start"
   mkdir -p $output/nanopore_ref_seq
   python $bin/python/extract_nanopore_ref.py $output/gmap_after_caun/uniq_gmap_after_caun.gff3 $output/CD-HIT-EST/lordec-correct_after_caun_cd_hit_est.fasta >$output/nanopore_ref_seq/nanopore_ref.fa
   ##transdecoder pep sequence
   /home/lfgu/tool/Trinity/trinityrnaseq_r20140413p1/trinity-plugins/transdecoder/TransDecoder -t $output/nanopore_ref_seq/nanopore_ref.fa  -S --reuse --workdir $output/nanopore_ref_seq
   mv nanopore_ref.fa.transdecoder.* $output/nanopore_ref_seq
   ##bowtie2 index
   bowtie2-build $output/nanopore_ref_seq/nanopore_ref.fa  nanopore_ref
   mv nanopore_ref.* $output/nanopore_ref_seq/
   #extract go term for biogo after blast2go
   python /home/lfgu/bigdata/hitseq/why/bin/python/blast2go_uniq.py $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.txt $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.cds $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.cds.uniq  $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.uniq $output/nanopore_ref_seq/Biological_Process.bingo $output/nanopore_ref_seq/Cellular_Component.bingo $output/nanopore_ref_seq/Molecular_Function.bingo > $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt
   #gmap cds_uniq
   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 60 -D $gmap/Cunninghamia_lanceolata/ --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f 2 -d Cunninghamia_lanceolata $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.cds.uniq >$output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.cds.uniq.gff3
   ##extract cds sequence from Genome
   python /home/lfgu/bigdata/hitseq/why/bin/python/extract_gff3_exon.py $Genome  $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.cds.uniq.gff3 > $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.cds.uniq.genome
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----nanopore seq-----'
fi


if [ ! -d "$output/cuffmerge_add_nanopore/RSEM" ]; then
        echo
        echo
        echo "[`date`] Beginning cuffmerge_add_nanopore/RSEM"
        echo '-----------------------------------------------'
        mkdir -p $output/cuffmerge_add_nanopore/RSEM
        mkdir -p $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/
        rsem-prepare-reference --bowtie2  $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas  $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore
        wait

        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Control_R1_rep1 $Control_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Control_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Control_R1_rep2 $Control_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Control_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Control_R1_rep3 $Control_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Control_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_1_R1_rep1 $Cold_1_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_1_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_1_R1_rep2 $Cold_1_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_1_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_1_R1_rep3 $Cold_1_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_1_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_2_R1_rep1 $Cold_2_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_2_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_2_R1_rep2 $Cold_2_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_2_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_2_R1_rep3 $Cold_2_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_2_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_3_R1_rep1 $Cold_3_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_3_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_3_R1_rep2 $Cold_3_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_3_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $Cold_3_R1_rep3 $Cold_3_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/Cold_3_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_1_R1_rep1 $PEG_1_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_1_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_1_R1_rep2 $PEG_1_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_1_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_1_R1_rep3 $PEG_1_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_1_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_2_R1_rep1 $PEG_2_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_2_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_2_R1_rep2 $PEG_2_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_2_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_2_R1_rep3 $PEG_2_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_2_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_3_R1_rep1 $PEG_3_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_3_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_3_R1_rep2 $PEG_3_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_3_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $PEG_3_R1_rep3 $PEG_3_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/PEG_3_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_1_R1_rep1 $SA_1_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_1_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_1_R1_rep2 $SA_1_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_1_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_1_R1_rep3 $SA_1_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_1_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_2_R1_rep1 $SA_2_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_2_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_2_R1_rep2 $SA_2_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_2_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_2_R1_rep3 $SA_2_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_2_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_3_R1_rep1 $SA_3_R2_rep1 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_3_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_3_R1_rep2 $SA_3_R2_rep2 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_3_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2-mismatch-rate 0.15 --bowtie2 $SA_3_R1_rep3 $SA_3_R2_rep3 $output/cuffmerge_add_nanopore/RSEM/index_real_nanopore_Seq/nanopore $output/cuffmerge_add_nanopore/RSEM/SA_3_rep3 &
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/cuffmerge_add_nanopore/RSEM/FPKM_plot" ]; then
        echo
        echo
        echo "[`date`] Beginning cuffmerge_add_nanopore/RSEM/FPKM_plot"
        echo '-----------------------------------------------'
        mkdir -p $output/cuffmerge_add_nanopore/RSEM/FPKM_plot
        for i in Control_rep1 Control_rep2 Control_rep3 Cold_1_rep1 Cold_1_rep2 Cold_1_rep3 Cold_2_rep1 Cold_2_rep2 Cold_2_rep3 Cold_3_rep1 Cold_3_rep2 Cold_3_rep3 PEG_1_rep1 PEG_1_rep2 PEG_1_rep3 PEG_2_rep1 PEG_2_rep2 PEG_2_rep3 PEG_3_rep1 PEG_3_rep2 PEG_3_rep3 SA_1_rep1 SA_1_rep2 SA_1_rep3 SA_2_rep1 SA_2_rep2 SA_2_rep3 SA_3_rep1 SA_3_rep2 SA_3_rep3
        do
        perl /home/lfgu/bigdata/hitseq/why/bin/tool/gra/gra/bin/RPKMplotV1.pl $output/cuffmerge_add_nanopore/RSEM/"$i".isoforms.results 7 $i >$output/cuffmerge_add_nanopore/RSEM/FPKM_plot/"$i".plot 2>$output/cuffmerge_add_nanopore/RSEM/FPKM_plot/"$i".stat
        done
        cat $output/cuffmerge_add_nanopore/RSEM/FPKM_plot/*.plot >$output/cuffmerge_add_nanopore/RSEM/FPKM_plot/all.plot
        Rscript /home/lfgu/tool/bin/hitseq/RPKMplotWithoutTitle.R $output/cuffmerge_add_nanopore/RSEM/FPKM_plot/fpkm.pdf $output/cuffmerge_add_nanopore/RSEM/FPKM_plot/all.plot
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/RSEM" ]; then
        echo
        echo
        echo "[`date`] Beginning RSEM"
        echo '-----------------------------------------------'
        mkdir -p $output/RSEM
        mkdir -p $output/RSEM/index/
        rsem-prepare-reference --bowtie2  $output/nanopore_ref_seq/nanopore_ref.fa  $output/RSEM/index/nanopore
        wait
     
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Control_R1_rep1 $Control_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/Control_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Control_R1_rep2 $Control_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/Control_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Control_R1_rep3 $Control_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/Control_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_1_R1_rep1 $Cold_1_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/Cold_1_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_1_R1_rep2 $Cold_1_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/Cold_1_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_1_R1_rep3 $Cold_1_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/Cold_1_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_2_R1_rep1 $Cold_2_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/Cold_2_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_2_R1_rep2 $Cold_2_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/Cold_2_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_2_R1_rep3 $Cold_2_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/Cold_2_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_3_R1_rep1 $Cold_3_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/Cold_3_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_3_R1_rep2 $Cold_3_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/Cold_3_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $Cold_3_R1_rep3 $Cold_3_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/Cold_3_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_1_R1_rep1 $PEG_1_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/PEG_1_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_1_R1_rep2 $PEG_1_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/PEG_1_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_1_R1_rep3 $PEG_1_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/PEG_1_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_2_R1_rep1 $PEG_2_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/PEG_2_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_2_R1_rep2 $PEG_2_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/PEG_2_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_2_R1_rep3 $PEG_2_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/PEG_2_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_3_R1_rep1 $PEG_3_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/PEG_3_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_3_R1_rep2 $PEG_3_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/PEG_3_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $PEG_3_R1_rep3 $PEG_3_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/PEG_3_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_1_R1_rep1 $SA_1_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/SA_1_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_1_R1_rep2 $SA_1_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/SA_1_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_1_R1_rep3 $SA_1_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/SA_1_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_2_R1_rep1 $SA_2_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/SA_2_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_2_R1_rep2 $SA_2_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/SA_2_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_2_R1_rep3 $SA_2_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/SA_2_rep3 &
        wait
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_3_R1_rep1 $SA_3_R2_rep1 $output/RSEM/index/nanopore $output/RSEM/SA_3_rep1 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_3_R1_rep2 $SA_3_R2_rep2 $output/RSEM/index/nanopore $output/RSEM/SA_3_rep2 &
        rsem-calculate-expression --paired-end -p 15 --bowtie2 $SA_3_R1_rep3 $SA_3_R2_rep3 $output/RSEM/index/nanopore $output/RSEM/SA_3_rep3 &
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

TRINITY_HOME=/home/lfgu/bigdata/hitseq/why/bin/tool/Trinityrnaseq-v2.6.6
if [ ! -d "$output/edgeR" ];then
   echo "[`date`]  start"
   echo '--------------------'
   mkdir -p $output/edgeR
   conda deactivate
   find ./output_nanopore/RSEM/ -maxdepth 2  -name "*.isoforms.results" | tee ./output_nanopore/RSEM/rsem.isoform.files
   $TRINITY_HOME/util/abundance_estimates_to_matrix.pl  --est_method RSEM  --quant_files $output/RSEM/rsem.isoform.files --out_prefix rsem  --gene_trans_map none  --name_sample_by_basedir
   #sort the columns
   python $bin/python/sort_RSEM_FPKM.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix  $output/RSEM/rsem.isoform.TMM.EXPR.sort.matrix  
   mv $output/rsem.isoform* $output/edgeR
   $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl -m $output/RSEM/rsem.isoform.counts.matrix -s $input/RNA_seq_samples.txt --method edgeR -o $output/edgeR --dispersion 0.1 --contrasts $input/contrasts.txt
   echo "[`date`]  finsh"
   echo '-----edgeR-----'
fi

TRINITY_HOME=/home/lfgu/bigdata/hitseq/why/bin/tool/Trinityrnaseq-v2.6.6
if [ ! -d "$output/cuffmerge_add_nanopore/edgeR" ];then
   echo "[`date`]  start"
   echo '--------------------'
   mkdir -p $output/cuffmerge_add_nanopore/edgeR
   conda deactivate
   find ./output_nanopore/cuffmerge_add_nanopore/RSEM/ -maxdepth 2  -name "*.isoforms.results" | tee ./output_nanopore/cuffmerge_add_nanopore/RSEM/rsem.isoform.files
   $TRINITY_HOME/util/abundance_estimates_to_matrix.pl  --est_method RSEM  --quant_files $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.files --out_prefix rsem  --gene_trans_map none
   mv rsem.isoform.* $output/cuffmerge_add_nanopore/RSEM/
   #sort the columns
   python $bin/python/sort_RSEM_FPKM.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix  $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix
   ##PCA
   more $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix |perl -ne 'chomp;if ($_=~/locus/){print "$_\n";}else{@l=split;next if (($l[1]<1)||($l[2]<1)||($l[3]<1)||($l[4]<1)||($l[5]<1)||($l[6]<1)||($l[7]<1)||($l[8]<1)||($l[9]<1)||($l[10]<1)||($l[11]<1)||($l[12]<1)||($l[13]<1)||($l[14]<1)||($l[15]<1)||($l[16]<1)||($l[17]<1)||($l[18]<1)||($l[19]<1)||($l[20]<1)||($l[21]<1)||($l[22]<1)||($l[23]<1)||($l[24]<1)||($l[25]<1)||($l[26]<1)||($l[27]<1)||($l[28]<1)||($l[29]<1)||($l[30]<1));print "$l[0]\t";for ($i=1;$i<$#l;$i++){$v=int($l[$i]);print "$v\t";};$v=int($l[30]);print "$v\n";}' > $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.PCA_hcluster.plot
   PCA_hcluster.R $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.PCA_hcluster.plot $output/FPKM/PCA_hcluster.pdf
   ##heatmap
   python $bin/python/SDX_nanopore_cuffmerge_FPKM_heatmap.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix > $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.heatmap
   Rscript $bin/R/heatmap_with_cluster.R $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.heatmap.pdf $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.heatmap
   Rscript $bin/R/heatmap_chengxu_without_cluster.R $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.heatmap.without.cluster.pdf $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix.heatmap
   ##
   mv $output/cuffmerge_add_nanopore/rsem.isoform* $output/cuffmerge_add_nanopore/edgeR
   $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl -m $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.counts.matrix -s $input/RNA_seq_samples.txt --method edgeR -o $output/cuffmerge_add_nanopore/edgeR --dispersion 0.1 --contrasts $input/contrasts.txt
   echo "[`date`]  finsh"
   echo '-----edgeR-----'
fi

if [ ! -d "$output/DEG" ];then
   echo "[`date`]  start"
   echo '--------------------'
   mkdir -p $output/DEG
   conda deactivate
   for i in Control
   do
   for j in Cold_1 Cold_2 Cold_3  PEG_1 PEG_2 PEG_3 SA_1 SA_2 SA_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt 
   done
   done

   
   for i in Cold_1
   do
   for j in Cold_2 Cold_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in Cold_2
   do
   for j in Cold_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in PEG_1
   do
   for j in PEG_2 PEG_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in PEG_2
   do
   for j in PEG_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in SA_1
   do
   for j in SA_2 SA_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt
   done
   done
   
   for i in SA_2
   do
   for j in SA_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/RSEM/rsem.isoform.TMM.EXPR.matrix $output/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt "$i" "$j" $output/DEG/"$i"_vs_"$j".txt
   done
   done
   wait

   ##split down and up regulate
   for i in Control
   do
   for j in Cold_1 Cold_2 Cold_3 PEG_1 PEG_2 PEG_3 SA_1 SA_2 SA_3
   do
   grep 'up'  $output/DEG/"$i"_vs_"$j".txt >  $output/DEG/"$i"_vs_"$j".up-regulation.xls
   grep 'down'  $output/DEG/"$i"_vs_"$j".txt >  $output/DEG/"$i"_vs_"$j".down-regulation.xls
   done
   done

   ##edgeR overlap
   python $bin/python/DEG_overlap.py $output/DEG/Control_vs_Cold_1.up-regulation.xls $output/DEG/Control_vs_Cold_2.up-regulation.xls $output/DEG/Control_vs_Cold_3.up-regulation.xls $output/DEG/Control_vs_PEG_1.up-regulation.xls $output/DEG/Control_vs_PEG_2.up-regulation.xls $output/DEG/Control_vs_PEG_3.up-regulation.xls $output/DEG/Control_vs_SA_1.up-regulation.xls $output/DEG/Control_vs_SA_2.up-regulation.xls $output/DEG/Control_vs_SA_3.up-regulation.xls $output/DEG/overlap/Cold_1_Cold_2_Cold_3_up_overlap.xls $output/DEG/overlap/PEG_1_PEG_2_PEG_3_up_overlap.xls $output/DEG/overlap/SA_1_SA_2_SA_3_up_overlap.xls $output/DEG/overlap/Cold_1_PEG_1_SA_1_up_overlap.xls $output/DEG/overlap/Cold_2_PEG_2_SA_2_up_overlap.xls $output/DEG/overlap/Cold_3_PEG_3_SA_3_up_overlap.xls $output/DEG/overlap/Cold_1_Cold_2_Cold_3_PEG_1_PEG_2_PEG_3_SA_1_SA_2_SA_3_up_overlap.xls
   python $bin/python/DEG_overlap.py $output/DEG/Control_vs_Cold_1.down-regulation.xls $output/DEG/Control_vs_Cold_2.down-regulation.xls $output/DEG/Control_vs_Cold_3.down-regulation.xls $output/DEG/Control_vs_PEG_1.down-regulation.xls $output/DEG/Control_vs_PEG_2.down-regulation.xls $output/DEG/Control_vs_PEG_3.down-regulation.xls $output/DEG/Control_vs_SA_1.down-regulation.xls $output/DEG/Control_vs_SA_2.down-regulation.xls $output/DEG/Control_vs_SA_3.down-regulation.xls $output/DEG/overlap/Cold_1_Cold_2_Cold_3_down_overlap.xls $output/DEG/overlap/PEG_1_PEG_2_PEG_3_down_overlap.xls $output/DEG/overlap/SA_1_SA_2_SA_3_down_overlap.xls $output/DEG/overlap/Cold_1_PEG_1_SA_1_down_overlap.xls $output/DEG/overlap/Cold_2_PEG_2_SA_2_down_overlap.xls $output/DEG/overlap/Cold_3_PEG_3_SA_3_down_overlap.xls $output/DEG/overlap/Cold_1_Cold_2_Cold_3_PEG_1_PEG_2_PEG_3_SA_1_SA_2_SA_3_down_overlap.xls

   ##add transcription factor
   python $bin/python/fir_SDX_stress_merge_all_DEG_add_TF.py Control Cold/PEG/SA 3 $output/DEG/ /home/lfgu/bigdata/hitseq/SDX_fir/all_TF/overlap_TF/nanopore/nanopore_overlap.txt ~/bigdata/Genome/chinaFir/Genome_v2/ath_orth/fir_cuffmerge_add_nanopore_ath_pep_30.txt >$output/DEG/merge_all_DEG.txt
   grep -v 'NotTF'  $output/DEG/merge_all_DEG.txt >$output/DEG/transcription_factor.txt
   echo "[`date`]  finsh"
   echo '-----DEG-----'
fi

if [ ! -d "$output/cuffmerge_add_nanopore/DEG" ];then
   echo "[`date`]  start"
   echo '--------------------'
   mkdir -p $output/cuffmerge_add_nanopore/DEG
   ##format Blast2Go
   #python /home/lfgu/bigdata/hitseq/why/bin/python/blast2go_uniq.py $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.cds $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.cds_uniq  $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.oneline_uniq $output/cuffmerge_add_nanopore/Biological_Process_bingo $output/cuffmerge_add_nanopore/Cellular_Component_bingo $output/cuffmerge_add_nanopore/Molecular_Function_bingo > $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt
   conda deactivate
   for i in Control
   do
   for j in Cold_1 Cold_2 Cold_3  PEG_1 PEG_2 PEG_3 SA_1 SA_2 SA_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix $output/cuffmerge_add_nanopore/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt "$i" "$j" $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt
   done
   done


   for i in Cold_1
   do
   for j in Cold_2 Cold_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix $output/cuffmerge_add_nanopore/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt "$i" "$j" $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in Cold_2
   do
   for j in Cold_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix $output/cuffmerge_add_nanopore/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt "$i" "$j" $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in PEG_1
   do
   for j in PEG_2 PEG_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix $output/cuffmerge_add_nanopore/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt "$i" "$j" $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in PEG_2
   do
   for j in PEG_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix $output/cuffmerge_add_nanopore/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt "$i" "$j" $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt
   done
   done

   for i in SA_1
   do
   for j in SA_2 SA_3
   do
   python $bin/python/run_trinity_DE_edgeR.py $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.matrix $output/cuffmerge_add_nanopore/edgeR/rsem.isoform.counts.matrix."$i"_vs_"$j".edgeR.DE_results $output/cuffmerge_add_nanopore/cuffmerge_add_nanopore_uniq.fas.transdecoder.pep.Blast2Go.uniq.txt "$i" "$j" $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt
   done
   done
   ##split down and up regulate
   for i in Control
   do
   for j in Cold_1 Cold_2 Cold_3 PEG_1 PEG_2 PEG_3 SA_1 SA_2 SA_3
   do
   grep 'up'  $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt >  $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".up-regulation.xls
   grep 'down'  $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".txt >  $output/cuffmerge_add_nanopore/DEG/"$i"_vs_"$j".down-regulation.xls
   done
   done

   ##edgeR overlap
   python $bin/python/DEG_overlap.py $output/cuffmerge_add_nanopore/DEG/Control_vs_Cold_1.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_Cold_2.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_Cold_3.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_PEG_1.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_PEG_2.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_PEG_3.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_SA_1.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_SA_2.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_SA_3.up-regulation.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_1_Cold_2_Cold_3_up_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/PEG_1_PEG_2_PEG_3_up_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/SA_1_SA_2_SA_3_up_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_1_PEG_1_SA_1_up_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_2_PEG_2_SA_2_up_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_3_PEG_3_SA_3_up_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_1_Cold_2_Cold_3_PEG_1_PEG_2_PEG_3_SA_1_SA_2_SA_3_up_overlap.xls
   python $bin/python/DEG_overlap.py $output/cuffmerge_add_nanopore/DEG/Control_vs_Cold_1.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_Cold_2.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_Cold_3.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_PEG_1.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_PEG_2.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_PEG_3.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_SA_1.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_SA_2.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/Control_vs_SA_3.down-regulation.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_1_Cold_2_Cold_3_down_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/PEG_1_PEG_2_PEG_3_down_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/SA_1_SA_2_SA_3_down_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_1_PEG_1_SA_1_down_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_2_PEG_2_SA_2_down_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_3_PEG_3_SA_3_down_overlap.xls $output/cuffmerge_add_nanopore/DEG/overlap/Cold_1_Cold_2_Cold_3_PEG_1_PEG_2_PEG_3_SA_1_SA_2_SA_3_down_overlap.xls

   ##add transcription factor
   #all DEG gene
   python $bin/python/fir_SDX_stress_merge_all_DEG_add_TF.py Control Cold/PEG/SA 3 $output/cuffmerge_add_nanopore/DEG/ /home/lfgu/bigdata/hitseq/SDX_fir/output_nanopore/cuffmerge_add_nanopore/TF/iTAK/tf_classification_onlyTF.txt ~/bigdata/Genome/chinaFir/Genome_v2/ath_orth/fir_cuffmerge_add_nanopore_ath_pep_30.txt  >$output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt
   #cold
   python $bin/python/fir_SDX_stress_merge_all_DEG_add_TF.py Control Cold 3 $output/cuffmerge_add_nanopore/DEG/ /home/lfgu/bigdata/hitseq/SDX_fir/output_nanopore/cuffmerge_add_nanopore/TF/iTAK/tf_classification_onlyTF.txt ~/bigdata/Genome/chinaFir/Genome_v2/ath_orth/fir_cuffmerge_add_nanopore_ath_pep_30.txt >$output/cuffmerge_add_nanopore/DEG/merge_cold_DEG.txt
   #PEG
   python $bin/python/fir_SDX_stress_merge_all_DEG_add_TF.py Control PEG 3 $output/cuffmerge_add_nanopore/DEG/ /home/lfgu/bigdata/hitseq/SDX_fir/output_nanopore/cuffmerge_add_nanopore/TF/iTAK/tf_classification_onlyTF.txt ~/bigdata/Genome/chinaFir/Genome_v2/ath_orth/fir_cuffmerge_add_nanopore_ath_pep_30.txt >$output/cuffmerge_add_nanopore/DEG/merge_PEG_DEG.txt
   #SA
   python $bin/python/fir_SDX_stress_merge_all_DEG_add_TF.py Control SA 3 $output/cuffmerge_add_nanopore/DEG/ /home/lfgu/bigdata/hitseq/SDX_fir/output_nanopore/cuffmerge_add_nanopore/TF/iTAK/tf_classification_onlyTF.txt ~/bigdata/Genome/chinaFir/Genome_v2/ath_orth/fir_cuffmerge_add_nanopore_ath_pep_30.txt >$output/cuffmerge_add_nanopore/DEG/merge_SA_DEG.txt
   grep -v 'NotTF'  $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt >$output/cuffmerge_add_nanopore/DEG/transcription_factor.txt
   ###DEG gene for DP_GP_cluster
   mkdir -p $output/cuffmerge_add_nanopore/DEG/DP_GP_cluster
   #all DEG gene
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_for_DP_GP_cluster.py $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix  >$output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/all/all_DEG_for_DP_GP_cluster.txt
   Rscript $bin/R/heatmap_with_cluster.R  $output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/all/all_DEG.heatmap.pdf $output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/all/all_DEG_for_DP_GP_cluster.txt
   #cold
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_for_DP_GP_cluster.py $output/cuffmerge_add_nanopore/DEG/merge_cold_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix  >$output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/cold/cold_DEG_for_DP_GP_cluster.txt
   #PEG
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_for_DP_GP_cluster.py $output/cuffmerge_add_nanopore/DEG/merge_PEG_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix  >$output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/PEG/PEG_DEG_for_DP_GP_cluster.txt
   #SA
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_for_DP_GP_cluster.py $output/cuffmerge_add_nanopore/DEG/merge_SA_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix  >$output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/SA/SA_DEG_for_DP_GP_cluster.txt
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate ageingml
   export PYTHONPATH=~/bigdata/hitseq/why/bin/tool/anaconda3/envs/ageingml/bin/
   python $bin/python/cor_miRNA_mRNA.py $output/RSEM/TCONS_00012497_ATAF1.txt $output/DEG/merge_all_DEG_bHLH_bZIP.txt > $output/DEG/ATAF1_bHLH_bZIP_cor.txt
   conda deactivate
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate ageingml
   export PYTHONPATH=~/bigdata/hitseq/why/bin/tool/anaconda3/envs/ageingml/bin/
   for i in all cold PEG SA
   do
   DP_GP_cluster.py -i $output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/"$i"/"$i"_DEG_for_DP_GP_cluster.txt -o $output/cuffmerge_add_nanopore/DEG/DP_GP_cluster/"$i"/"$i"_DEG -p pdf -n 20 --plot --fast &
   done
   conda deactivate
   #AP2 heatmap
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_Diff_TF.py $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix AP2 >$output/cuffmerge_add_nanopore/DEG/merge_all_DEG_AP2.txt
   Rscript $bin/R/heatmap_with_cluster.R $output/cuffmerge_add_nanopore/DEG/merge_all_DEG_AP2.pdf $output/cuffmerge_add_nanopore/DEG/merge_all_DEG_AP2.txt
   #NAC heatmap
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_Diff_TF.py $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix NAC >$output/cuffmerge_add_nanopore/DEG/merge_all_DEG_NAC.txt
   Rscript $bin/R/heatmap_with_cluster.R $output/cuffmerge_add_nanopore/DEG/merge_all_DEG_NAC.pdf $output/cuffmerge_add_nanopore/DEG/merge_all_DEG_NAC.txt
   #WRKY heatmap
   python $bin/python/chinese_fir_nanopore_add_cuffmerge_DEG_Diff_TF.py $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt $output/cuffmerge_add_nanopore/RSEM/rsem.isoform.TMM.EXPR.sort.matrix WRKY >$output/cuffmerge_add_nanopore/DEG/merge_all_DEG_WRKY.txt
   Rscript $bin/R/heatmap_with_cluster.R $output/cuffmerge_add_nanopore/DEG/merge_all_DEG_WRKY.pdf $output/cuffmerge_add_nanopore/DEG/merge_all_DEG_WRKY.txt
   ##down TF and up TF
   python $bin/python/SDX_cuffmerge_add_nanopore_DEG_bar.py $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.txt $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.bar.down.txt $output/cuffmerge_add_nanopore/DEG/merge_all_DEG.bar.up.txt 
   echo "[`date`]  finsh"
   echo '-----DEG-----'
fi
exit
if [ ! -d "$output/bingo" ];then
   echo "[`date`]  start"
   echo '--------------------'
   mkdir -p $output/bingo
   for i in Control_vs_Cold_1 Control_vs_Cold_2 Control_vs_Cold_3 Control_vs_PEG_1 Control_vs_PEG_2 Control_vs_PEG_3 Control_vs_SA_1 Control_vs_SA_2 Control_vs_SA_3 Cold_1_vs_Cold_2 Cold_2_vs_Cold_3 PEG_1_vs_PEG_2 PEG_2_vs_PEG_3 SA_1_vs_SA_2 SA_2_vs_SA_3
   do
       sed -n '2,$p'  $output/DEG/"$i".txt|cut -f 1 >$output/bingo/"$i".txt
   done
   echo "[`date`]  finsh"
   echo '-----bingo-----'
fi

output_small_RNA=/home/lfgu/bigdata/hitseq/SDX_fir/output_nanopore_small_RNA
output_degradome=/home/lfgu/bigdata/hitseq/SDX_fir/output_degradome/nanopore_degradome
if [ ! -d "$output_small_RNA/PCC" ];then
   echo "[`date`]  start"
   echo '--------------------'
   mkdir -p $output_small_RNA/PCC
   conda deactivate
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate ageingml
   export PYTHONPATH=~/bigdata/hitseq/why/bin/tool/anaconda3/envs/ageingml/bin/  
   python $bin/python/cor_miRNA_mRNA.py $output/RSEM/rsem.isoform.TMM.EXPR.sort.matrix $output_small_RNA/miRNA_RPM/merge_miRNA_rpm_repeat.txt >$output_small_RNA/PCC/cor_miRNA_mRNA.txt

   conda deactivate
   echo "[`date`]  finsh"
   echo '-----Pearson Correlation Coefficient-----'
fi

if [ ! -f ./output_degradome/nanopore_degradome/degradome_DEG_all_miR.xls ]; then
        echo "[`date`] degradome_DEG"
        echo '-----------------------------------------------'
        python $bin/python/cleaveland_merge_DEG.py $output/DEG/Control_vs_Cold_1.txt Cold_1 $output/DEG/Control_vs_Cold_2.txt Cold_2 $output/DEG/Control_vs_Cold_3.txt Cold_3 $output/DEG/Control_vs_PEG_1.txt PEG_1 $output/DEG/Control_vs_PEG_2.txt PEG_2 $output/DEG/Control_vs_PEG_3.txt PEG_3  $output/DEG/Control_vs_SA_1.txt SA_1 $output/DEG/Control_vs_SA_2.txt SA_2 $output/DEG/Control_vs_SA_3.txt SA_3 $output_small_RNA/edgeR/Control_vs_Cold_1.down.xls Cold_1_down $output_small_RNA/edgeR/Control_vs_Cold_1.up.xls Cold_1_up $output_small_RNA/edgeR/Control_vs_Cold_2.down.xls Cold_2_down $output_small_RNA/edgeR/Control_vs_Cold_2.up.xls Cold_2_up $output_small_RNA/edgeR/Control_vs_Cold_3.down.xls Cold_3_down $output_small_RNA/edgeR/Control_vs_Cold_3.up.xls Cold_3_up $output_small_RNA/edgeR/Control_vs_PEG_1.down.xls PEG_1_down $output_small_RNA/edgeR/Control_vs_PEG_1.up.xls PEG_1_up $output_small_RNA/edgeR/Control_vs_PEG_2.down.xls PEG_2_down $output_small_RNA/edgeR/Control_vs_PEG_2.up.xls PEG_2_up $output_small_RNA/edgeR/Control_vs_PEG_3.down.xls PEG_3_down $output_small_RNA/edgeR/Control_vs_PEG_3.up.xls PEG_3_up $output_small_RNA/edgeR/Control_vs_SA_1.down.xls SA_1_down $output_small_RNA/edgeR/Control_vs_SA_1.up.xls SA_1_up $output_small_RNA/edgeR/Control_vs_SA_2.down.xls SA_2_down $output_small_RNA/edgeR/Control_vs_SA_2.up.xls SA_2_up $output_small_RNA/edgeR/Control_vs_SA_3.down.xls SA_3_down $output_small_RNA/edgeR/Control_vs_SA_3.up.xls SA_3_up $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt ./output_degradome/nanopore_degradome/results_merge_all.txt $output_small_RNA/PCC/cor_miRNA_mRNA.txt  >> ./output_degradome/nanopore_degradome/degradome_DEG_all_miR.xls
        echo "[`date`] finish"
        echo '-----------------------------------------------'
fi

if [ ! -f ./output_degradome/nanopore_degradome/degradome_DEG.xls ]; then
        echo "[`date`] degradome_DEG"
        echo '-----------------------------------------------'
        python $bin/python/cleaveland_merge_DEG.py $output/DEG/Control_vs_Cold_1.txt Cold_1 $output/DEG/Control_vs_Cold_2.txt Cold_2 $output/DEG/Control_vs_Cold_3.txt Cold_3 $output/DEG/Control_vs_PEG_1.txt PEG_1 $output/DEG/Control_vs_PEG_2.txt PEG_2 $output/DEG/Control_vs_PEG_3.txt PEG_3  $output/DEG/Control_vs_SA_1.txt SA_1 $output/DEG/Control_vs_SA_2.txt SA_2 $output/DEG/Control_vs_SA_3.txt SA_3 $output/nanopore_ref_seq/nanopore_ref.fa.transdecoder.pep.blast2go.uniq.txt ./output_degradome/nanopore_degradome/results_merge.txt >> ./output_degradome/nanopore_degradome/degradome_DEG.xls
        echo "[`date`] finish"
        echo '-----------------------------------------------'
fi
exit
#cd-hit-est -i  $output/lordec-correct_after_caun/correct.fasta -o $output/CD-HIT-EST/lordec-correct_after_caun_cd_hit_est.fasta -c 0.9 -n 10 -d 0 -M 0 -T 0 -r 1
#cd_hit_est_clstr.pl $output/CD-HIT-EST/lordec-correct_after_caun_cd_hit_est.fasta.clstr >$output/CD-HIT-EST/lordec-correct_after_caun_cd_hit_est.fasta.clstr.format
exit

                cd-hit-est -i $output/lordec-correct/format_correct.fasta -o $output/CD-HIT-EST/fomat_correct_cd_hit_est.fasta -c 0.9 -n 10 -d 0 -M 0 -T 0 -r 1
                cd_hit_est_clstr.pl $output/CD-HIT-EST/fomat_correct_cd_hit_est.fasta.clstr >$output/CD-HIT-EST/fomat_correct_cd_hit_est.fasta.clstr.format

exit
~/bigdata/hitseq/why/bin/tool/ORFfinder -in $output/CD-HIT-EST/cd_hit_est.fasta -out $output/CD-HIT-EST/ORFfinder -ml 150
exit
if [ ! -d "$output/flash/" ]; then
        echo
        echo
        echo "[`date`] flash"
        echo '-----------------------------------------------'
        mkdir $output/flash
        source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate nanopore
        export PYTHONPATH=~/bigdata/hitseq/why/bin/tool/anaconda3/envs/nanopore/bin/
        flash  /home/lfgu/bigdata/hitseq/SDX_fir/input/RNA-seq-rep1/R1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/input/RNA-seq-rep1/R2.fastq -d $output/flash 
        echo "[`date`]  finsh"
        echo '-----lordec-correct-----'
fi

if [ ! -d "$output/proovread/" ]; then
        echo
        echo
        echo "[`date`] proovread"
        echo '-----------------------------------------------'
        mkdir $output/proovread
        /home/lfgu/bigdata/hitseq/why/bin/tool/proovread/bin/SeqChunker -s 20G -o $output/proovread/SeqChunk.%03d.fa $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fa
        echo "[`date`]  finsh"
        echo '-----SeqChunker-----'
fi
/home/lfgu/bigdata/hitseq/why/bin/tool/proovread/bin/proovread -l $output/proovread/SeqChunk.001.fa -s $output/flash/extendedFrags.fastq -s $output/flash/notCombined_1.fastq -s $output/flash/notCombined_2.fastq -t 50 --pre proovread-001 --overwrite --no-sampling 
exit
if [ ! -d "$output/lordec-correct" ]; then
   echo '---------------------------'
   echo "[`date`]  start"
   mkdir -p $output/lordec-correct
   echo '--------------------'
   #fastq_fasta_U2T
   python $bin/python/nanopore_fastq_fasta_U2T.py $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fq $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fa
   #lordec-correct
   lordec-correct -T 35 -k 21 -s 3 -S $output/lordec-correct/statistics.log -i $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fa -2 $input/all_RNA_seq.fastq -o $output/lordec-correct/correct.fasta &
   echo "[`date`]  finsh"
   echo '-----lordec-correct-----'
fi

if [ ! -d "$output/CD-HIT-EST/" ]; then
                echo
                echo
                echo "[`date`] CD-HIT-EST"
                echo '-----------------------------------------------'
                mkdir -p $output/CD-HIT-EST
                cd-hit-est -i $output/lordec-correct/correct.fasta -o $output/CD-HIT-EST/cd_hit_est.fasta -c 0.9 -n 10 -d 0 -M 20000 -T 0 -r 1
                cd_hit_est_clstr.pl $output/CD-HIT-EST/cd_hit_est.fasta.clstr >$output/CD-HIT-EST/cd_hit_est.fasta.clstr.format
                echo "[`date`] Run complete"
                echo '-----------------------------------------------'
fi
exit

if [ ! -d "$output/caun" ]; then
   echo '---------------------------'
   echo "[`date`]  start"
   mkdir -p $output/caun
   canu useGrid=false -correct gnuplotImageFormat=svg corOutCoverage=10000 mhapSensitivity=high corMinCoverage=0 correctedErrorRate=0.16 overlapper=minimap minReadLength=200 stopOnLowCoverage=0.01  corThreads=30 cormmapThreads=30 corovlThreads=30  minOverlapLength=100 genomeSize=10g -p cla_nanopore -d $output/caun -nanopore-raw $input/nanopore/Cla_SDX/Cla_SDX_mRNA/20190404_0918_MN29143_FAK54561_0668a83c/all_fastq_pass.fq batMemory=250 batThreads=40
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----caun correct-----'
fi

if [ ! -d "$output/porechop" ]; then
   echo '---------------------------'
   echo "[`date`]  start"
   mkdir -p $output/porechop
   porechop --format fasta -t 47 -i $output/caun/cla_nanopore.correctedReads.fasta -o $output/porechop/cla.choped.fasta
   echo '---------------------------'
   echo "[`date`]  finsh"
   echo '-----------porechop--------'
fi

if [ ! -d "$output/statistics" ]; then
   echo '---------------------------'
   echo '--------statistics---------'
   mkdir -p $output/statistics
   python $bin/python/pab_nanopore_statistics.py $input/nanopore/cla_nanopore.fa >$output/statistics/read_statistics.txt
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----statistics-----'
fi
exit
if [ ! -d "$output/lordec-correct_after_caun" ]; then
   echo '--------------------'
   echo '-----lordec-correct after caun-----'
   echo "[`date`]  start"
   mkdir -p $output/lordec-correct_after_caun
   cat /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Cold_1_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Cold_2_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Cold_3_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/Control_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/PEG_1_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/PEG_2_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/PEG_3_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/SA_1_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/SA_2_rep1.fastq /home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/SA_3_rep1.fastq >/home/lfgu/bigdata/hitseq/SDX_fir/output/fastq_ht2_filter/lordec-correct.fastq
   #T=threads k=kmer_len(big genome=21,small genome=19) s=solid_threshold(2 or 3)  S=out statistics file  i=long read FASTA/Q file  2=short read FASTA/Q file(s)  o=corrected_read_file
   lordec-correct -T 35 -k 21 -s 3 -S $output/lordec-correct_after_caun/statistics.log -i $output/porechop/cla.choped.fasta -2 $output/fastq_ht2_filter/lordec-correct.fastq -o $output/lordec-correct_after_caun/correct.fasta &
   wait
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----lordec-correct after caun-----'
fi

if [ ! -d "$output/gmap_after_caun" ]; then
   echo '-----------------------'
   echo '---------gmap after_caun----------'
   echo "[`date`]  start"
   mkdir -p $output/gmap_after_caun
   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -t 45 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab --no-chimeras  --cross-species  -n 1 --min-intronlength 30 --max-intronlength-middle 200000  -f samse -d pab $output/lordec-correct_after_caun/correct.fasta> $output/gmap_after_caun/pab.sam
   wait
   
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate python27
   export PYTHONPATH=/home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/envs/python27/bin/python 
   samtools view -bT $Genome -o $output/gmap_after_caun/pab.bam -@ 30 $output/gmap_after_caun/pab.sam
   samtools sort $output/gmap_after_caun/pab.bam $output/gmap_after_caun/pab.sort
   samtools index $output/gmap_after_caun/pab.sort
   /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/deactivate python27
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----gmap_after_caun-----'
fi

#if [ ! -d "$output/cupcake" ]; then
   echo '-----------------------'
   echo '---------cupcake----------'
   echo "[`date`]  start"
   /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/deactivate nanopore
  # mkdir $output/cupcake
   source /home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/bin/activate python27
   export PYTHONPATH=/home/lfgu/bigdata/hitseq/why/bin/tool/anaconda3/envs/python27/bin/python
   collapse_isoforms_by_sam.py --input $output/lordec-correct_after_caun/correct.fasta -s $output/gmap_after_caun/pab.sort.sam --dun-merge-5-shorter -o $output/cupcake/pab_cupcake 
   #filter_by_count.py $pref.collapsed --min_count=2 >filter_by_count.stdout
   #filter_away_subset.py $pref.collapsed >filter_away_subset.stdout
   #filter_away_subset.py $pref.collapsed.min_fl_2
   #cDNA_Cupcake for assembly evaluation using 5-hour timepoint
   #collapse_isoforms_by_sam.py -c 0.95 -i 0.95 --input $read1 -s $sortedsam --dun-merge-5-shorter -o $pref
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----cupcake-----'
#fi
exit
if [ ! -d "$output/lordec-correct" ]; then
   echo '--------------------'
   echo '-----lordec-correct-----'
   echo "[`date`]  start"
   mkdir -p $output/lordec-correct
   #T=threads k=kmer_len(big genome=21,small genome=19) s=solid_threshold(2 or 3)  S=out statistics file  i=long read FASTA/Q file  2=short read FASTA/Q file(s)  o=corrected_read_file
   lordec-correct -T 35 -k 21 -s 3 -S $output/lordec-correct/statistics.log -i $input/nanopore/pab_nanopore.fa -2 /home/lfgu/bigdata/hitseq/spruce_CircRNA/output_RNA_seq_SMART/fastq/Pab.fq -o $output/lordec-correct/correct.fasta &
   wait

   lordec-trim -i $output/lordec-correct/correct.fasta -o $output/lordec-correct/correct_trim.fasta &
   wait

   lordec-stat -T 35 -a 5 -i $input/nanopore/pab_nanopore.fa -2 /home/lfgu/bigdata/hitseq/spruce_CircRNA/output_RNA_seq_SMART/fastq/Pab.fq -k 21 -s 3 -S $output/lordec-correct/lordec_stat.txt
   wait
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----lordec-correct-----'
fi

if [ ! -d "$output/FLASh" ]; then
   echo '--------------------'
   echo '-----FLASh-----'
   echo "[`date`]  start"
   mkdir -p $output/FLASh 
   /home/lfgu/bigdata/hitseq/why/bin/tool/FLASH-1.2.11-Linux-x86_64/flash -M 100 -m 10 -t 40 /home/lfgu/bigdata/hitseq/spruce_CircRNA/output_RNA_seq_SMART/fastq/R1.fq /home/lfgu/bigdata/hitseq/spruce_CircRNA/output_RNA_seq_SMART/fastq/R2.fq -d $output/FLASh
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '-----lordec-correct-----'
fi


if [ ! -d "$output/gmap" ]; then
   echo '-----------------------'
   echo '---------gmap----------'
   echo "[`date`]  start"
   mkdir -p $output/gmap
   ###mapping all reads from nanopore sequence(without correct)
   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab -d pab  -t 60  --no-chimeras --cross-species  -n 1 -f 2 --min-intronlength 30 --max-intronlength-middle 200000 /home/lfgu/bigdata/hitseq/spruce_CircRNA/input/nanopore/pab_nanopore.fa > $output/gmap/pab_nanopore_uncorrect.gff3
   wait
  
   ###mapping trim reads from nanopore sequence after correct 
   #correct_trim.fasta is too large to use gmap
   split -l 2400000 $output/lordec-correct/correct_trim.fasta $output/lordec-correct/correct_trim
   wait
  
   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab -d pab  -t 60  --no-chimeras --cross-species  -n 1 -f 2 --min-intronlength 30 --max-intronlength-middle 200000 $output/lordec-correct/correct_trimaa > $output/gmap/pab_nanopore_correct_trim_aa.gff3
   wait

   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab -d pab  -t 60  --no-chimeras --cross-species  -n 1 -f 2 --min-intronlength 30 --max-intronlength-middle 200000 $output/lordec-correct/correct_trimab > $output/gmap/pab_nanopore_correct_trim_ab.gff3
   wait

   cat  $output/gmap/pab_nanopore_correct_trim_aa.gff3 $output/gmap/pab_nanopore_correct_trim_ab.gff3 >>$output/gmap/pab_nanopore_correct_trim.gff3
   wait
  
   rm $output/gmap/pab_nanopore_correct_trim_aa.gff3 $output/gmap/pab_nanopore_correct_trim_ab.gff3
   wait

   ###mapping reads from nanopore sequence after correct without trim
   #correct.fasta is too large to use gmap
   split -l 2000000 $output/lordec-correct/correct.fasta $output/lordec-correct/correct
   wait

   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab -d pab  -t 60  --no-chimeras --cross-species  -n 1 -f 2 --min-intronlength 30 --max-intronlength-middle 200000 $output/lordec-correct/correctaa >$output/gmap/pab_nanopore_correct_aa.gff3
   wait

   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab -d pab  -t 60  --no-chimeras --cross-species  -n 1 -f 2 --min-intronlength 30 --max-intronlength-middle 200000 $output/lordec-correct/correctab >$output/gmap/pab_nanopore_correct_ab.gff3
   wait

   /home/lfgu/bigdata/hitseq/why/bin/tool/gmap-2018-07-04/bin/gmapl.sse42 -D /home/lfgu/bigdata/Genome/Tree/spruce/gmap/pab -d pab  -t 60  --no-chimeras --cross-species  -n 1 -f 2 --min-intronlength 30 --max-intronlength-middle 200000 $output/lordec-correct/correctac >$output/gmap/pab_nanopore_correct_ac.gff3
   wait

   cat $output/gmap/pab_nanopore_correct_aa.gff3 $output/gmap/pab_nanopore_correct_ab.gff3 $output/gmap/pab_nanopore_correct_ac.gff3 >>$output/gmap/pab_nanopore_correct.gff3
   wait

   rm $output/gmap/pab_nanopore_correct_aa.gff3 $output/gmap/pab_nanopore_correct_ab.gff3 $output/gmap/pab_nanopore_correct_ac.gff3
   wait
   echo '--------------------'
   echo "[`date`]  finsh"
   echo '--------gmap--------'
fi



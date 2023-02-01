#!/usr/bin/bash

#To use absolute path for fastq files, if the folder was more than 1 sample, create a list as follow:
JOBID="Sample-100"
SAMPLE_IDS="Sample-100-GGTCAATA,Sample-100-ATCTTTAG,Sample-100-CAGAGGCC,Sample-100-TCAGCCGT"
TRANSCRIPTOME="/dfs3/swaruplab/smorabit/pipelines/sn-rna-seq/GRCh38.p12.mrna"
FASTQS="/dfs3/swaruplab/smorabit/data/AD_NucSeq_2019/SWA_13952_B01_NAN_Lane"

cellranger count --id=$JOBID \
                 --transcriptome=$TRANSCRIPTOME \
                 --fastqs=$FASTQS \
                 --sample=$SAMPLE_IDS \
                 --localcores=16 \
                 --localmem=64 \
                 --expect-cells=10000


echo "All processes for all samples were done !!"

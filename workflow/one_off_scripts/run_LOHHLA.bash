DIR_lohhla="/groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/run_lohhla"
patientId="ID13-Baseline-cfDNA-2017Mar03"
outputDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla/${patientId}"
BAMDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}"
normalBAMfile="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}/ID13-WBC.bam"
hlaPath="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/ID13-WBC/winners.hla.LOHHLA.txt"
hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta/ID13-WBC/hla_fasta.fa"
CopyNumLoc="/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza/${patientId}/${patientId}_alternative_solutions_LOHHLA.txt"

patientId="ID13-Evaluation-cfDNA-2017May23"
outputDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla/${patientId}"
BAMDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}"
normalBAMfile="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}/ID13-WBC.bam"
hlaPath="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/ID13-WBC/winners.hla.LOHHLA.txt"
hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta/ID13-WBC/hla_fasta.fa"
CopyNumLoc="/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza/${patientId}/${patientId}_alternative_solutions_LOHHLA.txt"

patientId="ID2-Baseline-cfDNA-2016Nov25"
outputDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla/${patientId}"
BAMDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}"
normalBAMfile="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}/ID2-WBC.bam"
hlaPath="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/ID2-WBC/winners.hla.LOHHLA.txt"
hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta/ID2-WBC/hla_fasta.fa"
CopyNumLoc="/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza/${patientId}/${patientId}_alternative_solutions_LOHHLA.txt"

patientId="ID2-Evaluation-cfDNA-2017Feb21"
outputDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla/${patientId}"
BAMDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}"
normalBAMfile="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}/ID2-WBC.bam"
hlaPath="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/ID2-WBC/winners.hla.LOHHLA.txt"
hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta/ID2-WBC/hla_fasta.fa"
CopyNumLoc="/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza/${patientId}/${patientId}_alternative_solutions_LOHHLA.txt"

patientId="T-001-1st-cfDNA"
outputDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla/${patientId}"
BAMDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}"
normalBAMfile="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}/T-001-WBC.bam"
hlaPath="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/T-001-WBC/winners.hla.LOHHLA.txt"
hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta/T-001-WBC/hla_fasta.fa"
CopyNumLoc="/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza/${patientId}/${patientId}_alternative_solutions_LOHHLA.txt"

patientId="T-001-3rd-cfDNA"
outputDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla/${patientId}"
BAMDir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}"
normalBAMfile="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/lohhla_bams/${patientId}/T-001-WBC.bam"
hlaPath="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/T-001-WBC/winners.hla.LOHHLA.txt"
hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta/T-001-WBC/hla_fasta.fa"
CopyNumLoc="/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza/${patientId}/${patientId}_alternative_solutions_LOHHLA.txt"


DIR_lohhla="/groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/run_lohhla"
Rscript /groups/wyattgrp/users/amunzur/gillian_proj/lohhla/LOHHLAscript.R \
--patientId ${patientId} \
--outputDir ${outputDir} \
--normalBAMfile ${normalBAMfile} \
--BAMDir ${BAMDir} \
--hlaPath ${hlaPath} \
--HLAfastaLoc ${hla_fasta} \
--CopyNumLoc ${CopyNumLoc} \
--mappingStep TRUE \
--minCoverageFilter 10 \
--fishingStep TRUE \
--cleanUp FALSE \
--gatkDir /home/amunzur/anaconda3/envs/polysolver/jar \
--novoDir /home/amunzur/anaconda3/envs/lohhla/bin  >> ${DIR_lohhla}/${patientId}.txt 2>&1 &

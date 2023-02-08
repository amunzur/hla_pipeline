ca polysolver
DIR_output="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types"
DIR_logs="/groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/call_hla_type"

find /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted -type f -name "*.bam" | grep -vE "(tmp|ID13-WBC|ID2-WBC|T-001-WBC)" | grep WBC | while read WBC_file; do
    WBC_name=$(basename ${WBC_file} .bam) # get the wbc name without the .bam extension
    mkdir ${DIR_output}/${WBC_name} # make the output dir for polysolver output
    rm ${DIR_logs}/${WBC_name}.txt # remove the old log files to avoid catting to the end of the file
    /home/amunzur/anaconda3/envs/polysolver/bin/shell_call_hla_type \
        ${WBC_file} \
        Unknown \
        1 \
        hg38 \
        STDFQ \
        0 \
        ${DIR_output}/${name}  >> ${DIR_logs}/${WBC_name}.txt 2>&1 &
done
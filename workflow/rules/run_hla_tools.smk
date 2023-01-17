# wildcard must be the WBC sample
rule call_hla_type:
    input:
        DIR_bams + "/subsetted_sorted/{wildcard}.bam"
    output:
        DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt", # the main output file needed
    params:
        DIR_output=DIR_results + "/polysolver/hla_types/{wildcard}" # output dir to feed into the command
    threads: 12
    conda: 
        "../envs/polysolver.yaml"
    shell:
            "shell_call_hla_type \
            {input} \
            Unknown \
            1 \
            hg38 \
            STDFQ \
            0 \
            {params}"

# wildcard must be the tumor name
rule call_hla_mutations_from_type:
    input:
        normal=lambda wildcards: get_WBC_bam("{wildcard}".format(wildcard=wildcards.wildcard)),
        tumor=DIR_bams + "/sorted/{wildcard}.bam",
        hla_results=DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt",
    output:
        directory(DIR_results + "/polysolver/hla_mutations/{wildcard}.bam"),
    conda: 
        "../envs/polysolver.yaml"
    threads: 12
    shell:
            "shell_call_hla_mutations_from_type \
            {input.normal} \
            {input.tumor} \
            {input.hla_results} \
            hg38 \
            STDFQ \
            {output}"


rule make_symlinks: 
    input:
        PATH_normalbam=lambda wildcards: get_WBC_bam("{wildcard}".format(wildcard=wildcards.wildcard)),
        PATH_tumorbam=DIR_bams + "/sorted/{wildcard}.bam",
        normalname=lambda wildcards: get_WBC_name("{wildcard}".format(wildcard=wildcards.wildcard)),
        tumorname=lambda wildcards: get_tumor_name("{wildcard}".format(wildcard=wildcards.wildcard)),
    output: 
        directory("results/data/bam/lohhla_bams/{wildcard}")
    run:
        shell(
            "ln -s {input.normalbam} "
        )


# rule run_lohhla:
#     input:
#         patientId=
#         normal=
#     output:
#         directory(DIR_results + "/lohhla/{wildcard}"),
#     conda:
#         "../envs/lohhla.yaml"
#     threads: 12
#     run:
#         shell(
#             'Rscript /groups/wyattgrp/users/amunzur/gillian_proj/lohhla/LOHHLAscript.R \
#             --patientId {input.patientId} \
#             --outputDir {output} \
#             --normalBAMfile {input.normal}} \
#             --BAMDir "${DIR_lohhla}/bams/${patientId}" \
#             --hlaPath "${DIR_lohhla}/hlas/${patientId}.txt" \
#             --HLAfastaLoc "${DIR_lohhla}/hla_fasta/${patientId}.fa" \
#             --CopyNumLoc "${DIR_lohhla}/solutions/${patientId}.txt" \
#             --mappingStep TRUE \
#             --minCoverageFilter 10 \
#             --fishingStep TRUE \
#             --cleanUp FALSE \
#             --gatkDir /home/amunzur/anaconda3/envs/polysolver/jar \
#             --novoDir /home/amunzur/anaconda3/envs/lohhla/bin >> ${DIR_lohhla}/logs/${patientId}.txt 2>&1 &'
#         )


# wildcard must be the WBC sample
rule call_hla_type:
    input:
        DIR_bams + "/sorted/{wildcard}.bam",
    output:
        DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt", # the main output file needed
    params:
        DIR_output=DIR_results + "/polysolver/hla_types/{wildcard}", # output dir to feed into the command
    threads: 12
    conda: 
        "../envs/polysolver.yaml"
    shell:
            "/home/amunzur/anaconda3/envs/polysolver/bin/shell_call_hla_type \
            {input} \
            Unknown \
            1 \
            hg38 \
            STDFQ \
            0 \
            {params}"

# wildcard is tumor bam
rule call_hla_mutations_from_type:
    input:
        NAME_tumor="{wildcard}",
        PATH_tumor=DIR_bams + "/sorted/{wildcard}.bam",
        PATH_hla_results=DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt",
    output:
        PATH_output=directory(DIR_results + "/polysolver/hla_mutations/{wildcard}.bam"),
    conda: 
        "../envs/polysolver.yaml"
    threads: 12
    shell:
        "shell_call_hla_mutations_from_type $(/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted/$(cat /groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_table.tsv | cut -f2,3 | grep {input.NAME_tumor} | cut -f1).bam) {input.PATH_tumor} {input.PATH_hla_results} hg38 STDFQ {output.PATH_output}"

# wildcard is WBC bam
rule generate_hla_fasta:
    input:
        polysolver_hla="/groups/wyattgrp/users/amunzur/gillian_proj/hla-polysolver/data/abc_complete.fasta",
        hla_types=DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt"
    output:
        hla_fasta=DIR_resources + "/lohhla_fasta/{wildcard}/hla_fasta.fa",
        tmp_hla_fasta=DIR_temp + "/{wildcard}/winner.hla.txt",
    run:
        shell('tr "\t" "\n" < <(cat {input.hla_types}  | cut -f 2,3)> {output.tmp_hla_fasta}')
        shell("grep --no-group-separator -A 1 -f {output.tmp_hla_fasta} {input.polysolver_hla} > {output.hla_fasta}")

rule index_hla_fasta:
    input:
        hla_fasta=DIR_resources+"/lohhla_fasta/{wildcard}/hla_fasta.fa"
    output:
        hla_fasta_idx=DIR_resources+"/lohhla_fasta/{wildcard}/hla_fasta.fa.fai"
    shell:
        "samtools faidx {input}"

# wildcard is tumor name
rule make_symlinks: 
    input:
        PATH_tumorbam=DIR_bams + "/sorted/{wildcard}.bam",
        PATH_tumorbai=DIR_bams + "/sorted/{wildcard}.bam.bai",
    output: 
        PATH_tumorbam=DIR_results + "/data/bam/lohhla_bams/{wildcard}/{wildcard}.bam",
        PATH_tumorbai=DIR_results + "/data/bam/lohhla_bams/{wildcard}/{wildcard}.bam.bai",
    run:
        import pandas as pd
        import os
        def make_symlinks(PATH_tumor):
            tumorbam = os.path.basename(PATH_tumor)
            paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_table.tsv", "\t", names = ["patient_id", "normalname", "tumorname", "wbcbam", "tumorbam"])
            paired_samples = paired_samples[["wbcbam", "tumorbam"]]

            mask = paired_samples["tumorbam"] == tumorbam
            wbcbam = paired_samples["wbcbam"][mask].tolist()[0]
            PATH_WBC = DIR_results + "/data/bam/sorted/" + wbcbam # wbc path
            DIR_output = DIR_results + "/data/bam/lohhla_bams/" + tumorbam.replace(".bam", "") #make the dir where symlinks for wbc and tumor bams will placed 

            # bams
            cmd1 = " ".join(["mkdir", DIR_output])
            cmd2 = " ".join(["ln -s", PATH_tumor, os.path.join(DIR_output, tumorbam)])
            cmd3 = " ".join(["ln -s", PATH_WBC, os.path.join(DIR_output, wbcbam)])
            # bais
            cmd4 = " ".join(["ln -s", PATH_tumor + ".bai", os.path.join(DIR_output, tumorbam + ".bai")])
            cmd5 = " ".join(["ln -s", PATH_WBC + ".bai", os.path.join(DIR_output, wbcbam + ".bai")])

            for cmd in [cmd1, cmd2, cmd3, cmd4, cmd5]:
                print(cmd)
                os.system(cmd)
        make_symlinks(input.PATH_tumorbam)

# rule modify_and_move_winners:
#     input:
#     output:
#     run:


#wildcard is wbc samples
# rule run_lohhla:
#     input:
#         tumorname="{wildcard}"
#         alternative_solutions=DIR_results + "/sequenza/{wildcard}/{wildcard}_alternative_solutions_LOHHLA.txt",
#         hla_fasta=DIR_resources + "/lohhla_fasta/{wildcard}/hla_fasta.fa",
#         lohhla_bam_dir=directory(DIR_results + "/data/bam/lohhla_bams/{wildcard}"),
#         hla_types= # winners from polysolvers
#     output:
#         directory(DIR_results + "/lohhla/{wildcard}"),
#     conda:
#         "../envs/lohhla.yaml"
#     threads: 12
#     run:
#         tumorbam = os.path.basename(PATH_tumor)
#         paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_table.tsv", "\t", names = ["patient_id", "normalname", "tumorname", "wbcbam", "tumorbam"])
#         paired_samples = paired_samples[["wbcbam", "tumorbam"]]

#         mask = paired_samples["tumorbam"] == tumorbam
#         wbcbam = paired_samples["wbcbam"][mask].tolist()[0]

#         cmd = 'Rscript /groups/wyattgrp/users/amunzur/gillian_proj/lohhla/LOHHLAscript.R --patientId input.tumorname --outputDir output --normalBAMfile input.normal --BAMDir input.lohhla_bam_dir --hlaPath input.hla_types --HLAfastaLoc input.hla_fasta --CopyNumLoc input.alternative_solutions --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir /home/amunzur/anaconda3/envs/polysolver/jar --novoDir /home/amunzur/anaconda3/envs/lohhla/bin'
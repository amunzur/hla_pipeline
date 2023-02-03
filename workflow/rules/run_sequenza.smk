# wildcard is tumor name
rule run_sequenza_utils:
    input:
        PATH_hg38=PATH_hg38,
        PATH_gcwig=PATH_gcwig,
        tumor_split_bams=return_split_chroms_tumor("{wildcard}"),
    params: 
        NAME_tumor="{wildcard}"
    output:
        seqz_chr1  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr1.seqz",
        seqz_chr2  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr2.seqz",
        seqz_chr3  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr3.seqz",
        seqz_chr4  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr4.seqz",
        seqz_chr5  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr5.seqz",
        seqz_chr6  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr6.seqz",
        seqz_chr7  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr7.seqz",
        seqz_chr8  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr8.seqz",
        seqz_chr9  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr9.seqz",
        seqz_chr10 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr10.seqz",
        seqz_chr11 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr11.seqz",
        seqz_chr12 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr12.seqz",
        seqz_chr13 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr13.seqz",
        seqz_chr14 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr14.seqz",
        seqz_chr15 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr15.seqz",
        seqz_chr16 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr16.seqz",
        seqz_chr17 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr17.seqz",
        seqz_chr18 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr18.seqz",
        seqz_chr19 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr19.seqz",
        seqz_chr20 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr20.seqz",
        seqz_chr21 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr21.seqz",
        seqz_chr22 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr22.seqz",
    conda:
        "../envs/sequenza.yaml"
    threads: 12
    shell:
        """
        DIR_bams_split='/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/split/'
        NAME_WBC=$(cat /groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_table.tsv | cut -f2,3 | grep {params.NAME_tumor} | cut -f1)
        PATH_WBC=$DIR_bams_split$NAME_WBC/
        
        DIR_results='/groups/wyattgrp/users/amunzur/hla_pipeline/results/sequenza'
        
        for CHROM in .REF_chr{{1..22}}.bam
        do
            sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n $PATH_WBC$NAME_WBC$CHROM -t $DIR_bams_split/{params.NAME_tumor}/{params.NAME_tumor}$CHROM -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > $DIR_results/{params.NAME_tumor}/{params.NAME_tumor}${{CHROM/.bam/.seqz}}
        done
        """

# puts the ploidy estimates file in such a format that lohhla would accept, minor reformatting
rule modify_ploidy_output:
    input:
        DIR_results + "/sequenza/{wildcard}/{wildcard}_alternative_solutions.txt",
    output:
        DIR_results + "/sequenza/{wildcard}/{wildcard}_alternative_solutions_LOHHLA.txt",
    params:
        samplename="{wildcard}"
    threads: 12
    run:
        shell('echo -e "Ploidy\ttumorPurity\ttumorPloidy\n{params.samplename}\t2\t$(tail -n +2 {input} | cut -f1 | head -n1)\t$(tail -n +2 {input} | cut -f2 | head -n1)" | awk -v OFS="\t" \'$1=$1\' > {output}')


# wildcard is tumor name
rule combine_seqz_files:
    input:
        return_seqz_per_chrom("{wildcard}", DIR_results),
    output:
        DIR_results + "/sequenza/{wildcard}/{wildcard}.seqz",
    shell: 
        "cat {input} > {output}"

# # wildcard is cfDNA samples, without the .bam extension
# rule run_sequenza_R:
# 	input:
# 		sample_name="{wildcard}",
# 		seqzfile=DIR_results + "/sequenza/{wildcard}/{wildcard}.seqz"
# 	output:
# 		DIR_output=directory(DIR_results + "/sequenza/{wildcard}")
# 		some_file=DIR_results + "/sequenza/{wildcard}" + "{wildcard}_alternative_solutions.txt" # just one of the sequenza output files is enough to show here
# 	conda:
# 		"../envs/sequenza.yaml"
# 	threads: 12
# 	shell:
# 		"Rscript /groups/wyattgrp/eritch/projects/ghent_m1rp/wxs/scripts/runseg/run_sequenza.R {input.seqzfile} {output.DIR_output} {input.sample_name}"

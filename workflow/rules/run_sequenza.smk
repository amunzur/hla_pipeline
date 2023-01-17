# wildcard is cfDNA name
rule run_sequenza_utils:
    input:
        wbc_chr1=lambda wildcards: expand("/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/split/{wildcard}/{wildcard}.REF_chr{{chr}}.bam".format(wildcard=wildcards.wildcard), chr=["1"]),
        PATH_hg38=PATH_hg38,
        PATH_gcwig=PATH_gcwig,
        tumor_chr1 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam",
        tumor_chr1_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam.bai",
    output:
        seqz_chr1  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr1.seqz",
    conda:
        "../envs/sequenza.yaml"
    threads: 12
    shell:
        'sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr1} -t {input.tumor_chr1} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output}'

lambda wildcards: expand("BAMs/{sample}_{{run}}.bam".format(sample=wildcards.sample), run=glob_wildcards("input_data/{sample}_{pair}_{{run}}.fastq.gz").format(sample=wildcards.sample, pair=PAIR1).run)





rule run_sequenza_utils:
	input:
		PATH_hg38=PATH_hg38,
		PATH_gcwig=PATH_gcwig,
		tumor_chr1  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam",
		tumor_chr2  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr2.bam",
		tumor_chr3  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr3.bam",
		tumor_chr4  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr4.bam",
		tumor_chr5  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr5.bam",
		tumor_chr6  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr6.bam",
		tumor_chr7  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr7.bam",
		tumor_chr8  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr8.bam",
		tumor_chr9  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr9.bam",
		tumor_chr10 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr10.bam",
		tumor_chr11 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr11.bam",
		tumor_chr12 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr12.bam",
		tumor_chr13 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr13.bam",
		tumor_chr14 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr14.bam",
		tumor_chr15 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr15.bam",
		tumor_chr16 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr16.bam",
		tumor_chr17 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr17.bam",
		tumor_chr18 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr18.bam",
		tumor_chr19 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr19.bam",
		tumor_chr20 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr20.bam",
		tumor_chr21 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr21.bam",
		tumor_chr22 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr22.bam",
		tumor_chrX  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrX.bam",
		tumor_chrY  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrY.bam",
		tumor_chr1_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam.bai",
		tumor_chr2_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr2.bam.bai",
		tumor_chr3_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr3.bam.bai",
		tumor_chr4_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr4.bam.bai",
		tumor_chr5_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr5.bam.bai",
		tumor_chr6_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr6.bam.bai",
		tumor_chr7_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr7.bam.bai",
		tumor_chr8_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr8.bam.bai",
		tumor_chr9_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr9.bam.bai",
		tumor_chr10_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr10.bam.bai",
		tumor_chr11_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr11.bam.bai",
		tumor_chr12_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr12.bam.bai",
		tumor_chr13_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr13.bam.bai",
		tumor_chr14_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr14.bam.bai",
		tumor_chr15_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr15.bam.bai",
		tumor_chr16_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr16.bam.bai",
		tumor_chr17_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr17.bam.bai",
		tumor_chr18_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr18.bam.bai",
		tumor_chr19_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr19.bam.bai",
		tumor_chr20_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr20.bam.bai",
		tumor_chr21_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr21.bam.bai",
		tumor_chr22_bai = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr22.bam.bai",
		tumor_chrX_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrX.bam.bai",
		tumor_chrY_bai  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrY.bam.bai",
		wbc_chr1  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "1", DIR_bams),
		wbc_chr2  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "2", DIR_bams),
		wbc_chr3  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "3", DIR_bams),
		wbc_chr4  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "4", DIR_bams),
		wbc_chr5  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "5", DIR_bams),
		wbc_chr6  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "6", DIR_bams),
		wbc_chr7  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "7", DIR_bams),
		wbc_chr8  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "8", DIR_bams),
		wbc_chr9  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "9", DIR_bams),
		wbc_chr10 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "10", DIR_bams),
		wbc_chr11 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "11", DIR_bams),
		wbc_chr12 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "12", DIR_bams),
		wbc_chr13 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "13", DIR_bams),
		wbc_chr14 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "14", DIR_bams),
		wbc_chr15 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "15", DIR_bams),
		wbc_chr16 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "16", DIR_bams),
		wbc_chr17 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "17", DIR_bams),
		wbc_chr18 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "18", DIR_bams),
		wbc_chr19 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "19", DIR_bams),
		wbc_chr20 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "20", DIR_bams),
		wbc_chr21 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "21", DIR_bams),
		wbc_chr22 = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "22", DIR_bams),
		wbc_chrX  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "X", DIR_bams),
		wbc_chrY  = lambda wildcards: get_WBC_chrom("{wildcard}".format(wildcard=wildcards.wildcard), "Y", DIR_bams),
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
		seqz_chrX  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chrX.seqz",
		seqz_chrY  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chrY.seqz",
	conda:
		"../envs/sequenza.yaml"
	threads: 12
	shell:
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr1} -t {input.tumor_chr1} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr1}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr2} -t {input.tumor_chr2} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr2}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr3} -t {input.tumor_chr3} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr3}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr4} -t {input.tumor_chr4} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr4}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr5} -t {input.tumor_chr5} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr5}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr6} -t {input.tumor_chr6} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr6}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr7} -t {input.tumor_chr7} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr7}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr8} -t {input.tumor_chr8} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr8}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr9} -t {input.tumor_chr9} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr9}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr10} -t {input.tumor_chr10} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr10}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr11} -t {input.tumor_chr11} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr11}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr12} -t {input.tumor_chr12} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr12}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr13} -t {input.tumor_chr13} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr13}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr14} -t {input.tumor_chr14} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr14}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr15} -t {input.tumor_chr15} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr15}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr16} -t {input.tumor_chr16} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr16}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr17} -t {input.tumor_chr17} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr17}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr18} -t {input.tumor_chr18} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr18}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr19} -t {input.tumor_chr19} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr19}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr20} -t {input.tumor_chr20} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr20}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr21} -t {input.tumor_chr21} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr21}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chr22} -t {input.tumor_chr22} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chr22}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chrX} -t {input.tumor_chrX} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chrX}")
		shell("sequenza-utils bam2seqz --fasta {input.PATH_hg38} -n {input.wbc_chrY} -t {input.tumor_chrY} -gc {input.PATH_gcwig} --het 0.4 -N 40 |  grep -v ^K | grep -v ^G | grep -v ^M | (sed -u 1q; sort -k1,1V -k2,2n) > {output.seqz_chrY}")


# wildcard is tumor name
# rule combine_seqz_files:
# 	input:
# 		seqz_chr1  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr1.seqz",
# 		seqz_chr2  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr2.seqz",
# 		seqz_chr3  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr3.seqz",
# 		seqz_chr4  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr4.seqz",
# 		seqz_chr5  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr5.seqz",
# 		seqz_chr6  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr6.seqz",
# 		seqz_chr7  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr7.seqz",
# 		seqz_chr8  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr8.seqz",
# 		seqz_chr9  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr9.seqz",
# 		seqz_chr10 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr10.seqz",
# 		seqz_chr11 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr11.seqz",
# 		seqz_chr12 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr12.seqz",
# 		seqz_chr13 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr13.seqz",
# 		seqz_chr14 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr14.seqz",
# 		seqz_chr15 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr15.seqz",
# 		seqz_chr16 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr16.seqz",
# 		seqz_chr17 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr17.seqz",
# 		seqz_chr18 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr18.seqz",
# 		seqz_chr19 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr19.seqz",
# 		seqz_chr20 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr20.seqz",
# 		seqz_chr21 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr21.seqz",
# 		seqz_chr22 = DIR_results + "/sequenza/{wildcard}/{wildcard}_chr22.seqz",
# 		seqz_chrX  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chrX.seqz",
# 		seqz_chrY  = DIR_results + "/sequenza/{wildcard}/{wildcard}_chrY.seqz",
# 	output:
# 		seqz_combined = DIR_results + "/sequenza/{wildcard}/{wildcard}.seqz",
# 	shell: 
# 		"cat {input} > {output}"
# # {input.seqz_chr1} {input.seqz_chr2} {input.seqz_chr3} {input.seqz_chr4} {input.seqz_chr5} {input.seqz_chr6} {input.seqz_chr7} {input.seqz_chr8} {input.seqz_chr9} {input.seqz_chr10} {input.seqz_chr11} {input.seqz_chr12} {input.seqz_chr13} {input.seqz_chr14} {input.seqz_chr15} {input.seqz_chr16} {input.seqz_chr17} {input.seqz_chr18} {input.seqz_chr19} {input.seqz_chr20} {input.seqz_chr21} {input.seqz_chr22} {input.seqz_chrX} {input.seqz_chrY}  > {output}

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

rule run_fastqc_raw:
    input:
        DIR_fastq + "/raw/{wildcard}.fq",
    output:
        output_zip=DIR_fastqc_raw + "/{wildcard}_fastqc.zip",
        output_html=DIR_fastqc_raw + "/{wildcard}_fastqc.html",
    threads: 5
    shell:
        "/home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output.output_zip}`"

# mask low quality bases in fastq files
rule mask_fastq:
    input:
        DIR_fastq + "/raw/{wildcard}.fq",
    output:
        DIR_fastq + "/masked/{wildcard}.fq",
    params:
        min_base_quality=20,
    run:
        shell(
            "/groups/wyattgrp/users/amunzur/gene_panel_pipeline/dependencies/fasta mask by quality {input} {params.min_base_quality} > {output}"
        )

rule trim_fastq:
    input:
        R1=DIR_fastq + "/masked/{wildcard}_1.fq",
        R2=DIR_fastq + "/masked/{wildcard}_2.fq",
    output:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq",
        html_report="results/reports/fastp/{wildcard}.html",
        json_report="results/reports/fastp/{wildcard}.json",
    threads: 12
    params:
        minimum_read_length=50,
    run:
        shell(
            "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} \
        --length_required {params.minimum_read_length} \
        --disable_trim_poly_g \
        --disable_quality_filtering \
        --html {output.html_report} \
        --json {output.json_report} \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --thread {threads}"
        )
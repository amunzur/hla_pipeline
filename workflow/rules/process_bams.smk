# Generate mapped bam
rule mapBAM:
    input:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/raw/{wildcard}.sam")
    threads: 12
    run:
        shell("bwa mem {params.PATH_hg38} {input.R1} {input.R2} -t {threads} > {output}")

rule samtobam:
    input:
        DIR_bams + "/raw/{wildcard}.sam"
    output:
        temp(DIR_bams + "/raw/{wildcard}.bam")
    run:
        shell("samtools view -b -h -o {output} {input}")

# Add read groups to the mapped bam file
rule addRG:
    input:
        DIR_bams + "/raw/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        sample="{wildcard}",
    output:
        temp(DIR_bams + "/rg/{wildcard}.bam")
    threads: 12
    run:
        shell(
            "picard AddOrReplaceReadGroups -Xmx20G I={input} O={output} RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1"
        )

rule fixmate:
    input:
        DIR_bams + "/rg/{wildcard}.bam",
    output:
        temp(DIR_bams + "/fixmate/{wildcard}.bam")
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input} O={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        )

rule markDuplicates:
    input: 
        DIR_bams + "/fixmate/{wildcard}.bam"
    output:
        bam=temp(DIR_bams + "/markDuplicates/{wildcard}.bam"),
        metrics=DIR_metrics + "/markDuplicates/{wildcard}.txt"
    run: 
        shell(
            "picard -Xmx40g MarkDuplicates  I={input} O={output.bam} M={output.metrics}"
        )

rule proper_pairs:
    input:
        DIR_bams + "/markDuplicates/{wildcard}.bam"
    output:
        temp(DIR_bams + "/proper_pair/{wildcard}.bam")
    threads: 12
    run:
        shell(
            "sambamba view {input} -F 'proper_pair' -t {threads} -f bam -l 0 -o {output}"
        )

rule sort_subseted_to_proper_pairs:
    input:
        DIR_bams + "/proper_pair/{wildcard}.bam"
    output:
        DIR_bams + "/sorted/{wildcard}.bam"
    threads: 12
    run:
        shell(
            "samtools sort {input} -o {output}"
        )

rule index_sorted:
    input:
        DIR_bams + "/sorted/{wildcard}.bam"
    output:
        DIR_bams + "/sorted/{wildcard}.bam.bai"
    run:
        shell(
            "samtools index {input}"
        )

# Remove unwanted contigs
rule filter_bam:
    input:
        bam=DIR_bams + "/sorted/{wildcard}.bam",
        bai=DIR_bams + "/sorted/{wildcard}.bam.bai"
    output:
        bam=DIR_bams + "/subsetted/{wildcard}.bam",
    run:
        shell(
            "samtools view {input.bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -bh > {output.bam}"
            )

rule sort_filtered:
    input:
        bam=DIR_bams + "/subsetted/{wildcard}.bam",
    output:
        bam=DIR_bams + "/subsetted_sorted/{wildcard}.bam",
    run:
        shell(
            "samtools sort {input} -o {output}"
            )

rule index_filtered:
    input:
        DIR_bams + "/subsetted_sorted/{wildcard}.bam",
    output:
        DIR_bams + "/subsetted_sorted/{wildcard}.bam.bai",
    run:
        shell(
            "samtools index {input}"
            )

rule split_bam:
    input:
        bam=DIR_bams + "/subsetted_sorted/{wildcard}.bam",
        bai=DIR_bams + "/subsetted_sorted/{wildcard}.bam.bai",
    output:
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr1.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr2.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr3.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr4.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr5.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr6.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr7.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr8.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr9.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr10.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr11.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr12.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr13.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr14.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr15.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr16.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr17.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr18.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr19.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr20.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr21.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr22.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chrX.bam",
        DIR_bams + "/subsetted_sorted/{wildcard}.REF_chrY.bam",
        # expand(DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr{chr_values}.bam", wildcard=samples, chr_values=["1", "2"]),
        # DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr{{X,Y}}.bam"
    shell:
        "bamtools split -in {input.bam} -reference"

rule move_split_bams:
    input:
        chr1  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr1.bam",
        chr2  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr2.bam",
        chr3  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr3.bam",
        chr4  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr4.bam",
        chr5  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr5.bam",
        chr6  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr6.bam",
        chr7  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr7.bam",
        chr8  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr8.bam",
        chr9  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr9.bam",
        chr10 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr10.bam",
        chr11 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr11.bam",
        chr12 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr12.bam",
        chr13 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr13.bam",
        chr14 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr14.bam",
        chr15 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr15.bam",
        chr16 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr16.bam",
        chr17 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr17.bam",
        chr18 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr18.bam",
        chr19 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr19.bam",
        chr20 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr20.bam",
        chr21 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr21.bam",
        chr22 = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chr22.bam",
        chrX  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chrX.bam",
        chrY  = DIR_bams + "/subsetted_sorted/{wildcard}.REF_chrY.bam",
    output:
        chr1  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam",
        chr2  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr2.bam",
        chr3  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr3.bam",
        chr4  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr4.bam",
        chr5  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr5.bam",
        chr6  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr6.bam",
        chr7  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr7.bam",
        chr8  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr8.bam",
        chr9  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr9.bam",
        chr10 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr10.bam",
        chr11 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr11.bam",
        chr12 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr12.bam",
        chr13 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr13.bam",
        chr14 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr14.bam",
        chr15 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr15.bam",
        chr16 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr16.bam",
        chr17 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr17.bam",
        chr18 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr18.bam",
        chr19 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr19.bam",
        chr20 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr20.bam",
        chr21 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr21.bam",
        chr22 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr22.bam",
        chrX  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrX.bam",
        chrY  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrY.bam",
    run:
        shell('mv {input.chr1} {output.chr1}')
        shell('mv {input.chr2} {output.chr2}')
        shell('mv {input.chr3} {output.chr3}')
        shell('mv {input.chr4} {output.chr4}')
        shell('mv {input.chr5} {output.chr5}')
        shell('mv {input.chr6} {output.chr6}')
        shell('mv {input.chr7} {output.chr7}')
        shell('mv {input.chr8} {output.chr8}')
        shell('mv {input.chr9} {output.chr9}')
        shell('mv {input.chr10} {output.chr10}')
        shell('mv {input.chr11} {output.chr11}')
        shell('mv {input.chr12} {output.chr12}')
        shell('mv {input.chr13} {output.chr13}')
        shell('mv {input.chr14} {output.chr14}')
        shell('mv {input.chr15} {output.chr15}')
        shell('mv {input.chr16} {output.chr16}')
        shell('mv {input.chr17} {output.chr17}')
        shell('mv {input.chr18} {output.chr18}')
        shell('mv {input.chr19} {output.chr19}')
        shell('mv {input.chr20} {output.chr20}')
        shell('mv {input.chr21} {output.chr21}')
        shell('mv {input.chr22} {output.chr22}')
        shell('mv {input.chrX} {output.chrX}')
        shell('mv {input.chrY} {output.chrY}')

rule index_split_bams:
    input:
        chr1  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam",
        chr2  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr2.bam",
        chr3  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr3.bam",
        chr4  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr4.bam",
        chr5  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr5.bam",
        chr6  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr6.bam",
        chr7  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr7.bam",
        chr8  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr8.bam",
        chr9  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr9.bam",
        chr10 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr10.bam",
        chr11 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr11.bam",
        chr12 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr12.bam",
        chr13 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr13.bam",
        chr14 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr14.bam",
        chr15 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr15.bam",
        chr16 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr16.bam",
        chr17 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr17.bam",
        chr18 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr18.bam",
        chr19 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr19.bam",
        chr20 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr20.bam",
        chr21 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr21.bam",
        chr22 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr22.bam",
        chrX  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrX.bam",
        chrY  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrY.bam",
    output:
        chr1  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr1.bam.bai",
        chr2  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr2.bam.bai",
        chr3  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr3.bam.bai",
        chr4  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr4.bam.bai",
        chr5  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr5.bam.bai",
        chr6  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr6.bam.bai",
        chr7  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr7.bam.bai",
        chr8  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr8.bam.bai",
        chr9  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr9.bam.bai",
        chr10 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr10.bam.bai",
        chr11 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr11.bam.bai",
        chr12 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr12.bam.bai",
        chr13 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr13.bam.bai",
        chr14 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr14.bam.bai",
        chr15 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr15.bam.bai",
        chr16 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr16.bam.bai",
        chr17 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr17.bam.bai",
        chr18 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr18.bam.bai",
        chr19 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr19.bam.bai",
        chr20 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr20.bam.bai",
        chr21 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr21.bam.bai",
        chr22 = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chr22.bam.bai",
        chrX  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrX.bam.bai",
        chrY  = DIR_bams + "/split/{wildcard}/{wildcard}.REF_chrY.bam.bai",
    run:
        shell('samtools index {input.chr1}')
        shell('samtools index {input.chr2}')
        shell('samtools index {input.chr3}')
        shell('samtools index {input.chr4}')
        shell('samtools index {input.chr5}')
        shell('samtools index {input.chr6}')
        shell('samtools index {input.chr7}')
        shell('samtools index {input.chr8}')
        shell('samtools index {input.chr9}')
        shell('samtools index {input.chr10}')
        shell('samtools index {input.chr11}')
        shell('samtools index {input.chr12}')
        shell('samtools index {input.chr13}')
        shell('samtools index {input.chr14}')
        shell('samtools index {input.chr15}')
        shell('samtools index {input.chr16}')
        shell('samtools index {input.chr17}')
        shell('samtools index {input.chr18}')
        shell('samtools index {input.chr19}')
        shell('samtools index {input.chr20}')
        shell('samtools index {input.chr21}')
        shell('samtools index {input.chr22}')
        shell('samtools index {input.chrX}')
        shell('samtools index {input.chrY}')


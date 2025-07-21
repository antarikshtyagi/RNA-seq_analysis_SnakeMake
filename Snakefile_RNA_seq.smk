SAMPLES, = glob_wildcards("raw_fastq/{sample}_R1.fastq.gz")

# path to track and reference
GTF = 'mm10_files/genes.gtf'
GENOME = 'mm10_files/genome.fa'
HISAT_INDEX = 'mm10_files/Hisat2Index/genome'

rule all:
    input:
        expand(["raw_fastqc/{sample}_R1_fastqc.html", "raw_fastqc/{sample}_R2_fastqc.html"], sample = SAMPLES),
        expand(["trimmed_fastqc/{sample}_paired_R1_fastqc.html", "trimmed_fastqc/{sample}_paired_R2_fastqc.html"], sample = SAMPLES),
        expand("bigwig/{sample}.bw", sample=SAMPLES),
        "RNA-Seq/Counts/featureCounts_gene_counts.txt"

rule raw_fastqc:
    input:
        expand(["raw_fastq/{sample}_R1.fastq.gz", "raw_fastq/{sample}_R2.fastq.gz"], sample = SAMPLES)
    output:
        expand(["raw_fastqc/{sample}_R1_fastqc.html", "raw_fastqc/{sample}_R2_fastqc.html"], sample = SAMPLES)
    resources: time_min=2880, mem_mb=20000, cpus=1
    params:
        out = "raw_fastqc",
        partition = "talon"
    shell:
        "fastqc -o {params.out} -t {resources.cpus} {input}"


rule trim_adapt:
    input:
        fq1 = ["raw_fastq/{sample}_R1.fastq.gz"],
        fq2 = ["raw_fastq/{sample}_R2.fastq.gz"]
    output:
        tr1 = ["trimmed_fastq/{sample}_paired_R1.fastq.gz"],
        tr2 = ["trimmed_fastq/{sample}_paired_R2.fastq.gz"],
        un1 = ["unpaired_fastq/{sample}_unpaired_R1.fastq.gz"],
        un2 = ["unpaired_fastq/{sample}_unpaired_R2.fastq.gz"],
        st = ["trimmed_fastq/{sample}_trimStats.txt"]
    resources: time_min=2880, mem_mb=20000, cpus=1
    params: partition = "talon"
    shell:
        "trimmomatic PE -threads {resources.cpus} -phred33 "
        "{input.fq1} {input.fq2} {output.tr1} {output.un1} {output.tr2} {output.un2} "
        "ILLUMINACLIP:./trimmomatic_adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:20 MINLEN:12 2> {output.st}"

rule trimmed_fastqc:
    input:
        expand(["trimmed_fastq/{sample}_paired_R1.fastq.gz", "trimmed_fastq/{sample}_paired_R2.fastq.gz"], sample = SAMPLES)
    output:
        expand(["trimmed_fastqc/{sample}_paired_R1_fastqc.html", "trimmed_fastqc/{sample}_paired_R2_fastqc.html"], sample = SAMPLES)
    resources: time_min=2880, mem_mb=20000, cpus=1
    params:
        out = "trimmed_fastqc",
        partition = "talon"
    shell:
        "fastqc -o {params.out} -t {resources.cpus} {input}"

rule hisat_map:
    input:
        fq1 = ["trimmed_fastq/{sample}_paired_R1.fastq.gz"],
        fq2 = ["trimmed_fastq/{sample}_paired_R2.fastq.gz"]
    output:
        bam = "hisat_mapped/{sample}.bam",
        STATS = "hisat_mapped/{sample}_mapping_stats.txt"
    resources: time_min=2880, mem_mb=20000, cpus=4
    params: partition = "talon"
    shell:
        "hisat2 -p {resources.cpus} -x {HISAT_INDEX} -1 {input.fq1} -2 {input.fq2} --new-summary --summary-file {output.STATS} --new-summary | "
        "samtools view -@ {resources.cpus} -S -b - > {output.bam}"

rule sam_sort:
    input:
        rules.hisat_map.output.bam
    output:
        "sorted_bam/{sample}_sorted.bam"
    resources: time_min=2880, mem_mb=20000, cpus=2
    params: partition = "talon"
    shell:
        "samtools sort -@ {resources.cpus} {input} -o {output}"

rule samtools_view_pair:
    input:
        rules.sam_sort.output
    output:
        "samtools_view_pair/{sample}_properly_paired_sorted.bam"
    resources: time_min=2880, mem_mb=20000, cpus=2
    params: partition = "talon"
    shell:
        "samtools view -@ {resources.cpus} -b -f 0x2 {input} > {output}"

rule sam_index:
    input:
        rules.samtools_view_pair.output
    output:
        "samtools_view_pair/{sample}_properly_paired_sorted.bam.bai"
    resources: time_min=2880, mem_mb=20000, cpus=1
    params: partition = "talon"
    shell:
        "samtools index -@ {resources.cpus} {input} {output}"

rule bamcoverage:
    input:
        bam = rules.samtools_view_pair.output,
        bai = rules.sam_index.output
    output:
        "bigwig/{sample}.bw"
    resources: time_min=2880, mem_mb=20000, cpus=4
    params: partition = "talon"
    shell:
        "bamCoverage -p {resources.cpus} -b {input.bam} --normalizeUsing CPM -of bigwig -o {output}"

rule featurecounts:
    input:
        bams = expand("samtools_view_pair/{sample}_properly_paired_sorted.bam", sample=SAMPLES)
    output:
        "RNA-Seq/Counts/featureCounts_gene_counts.txt"
    resources: time_min=2880, mem_mb=20000, cpus=24
    params: partition = "talon"
    shell:
        "featureCounts -p -T {resources.cpus} -a {GTF} -o {output} {input.bams}"

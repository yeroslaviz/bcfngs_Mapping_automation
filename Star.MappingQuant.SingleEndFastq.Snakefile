# the mpping file for either paired-end or songle-end mapping step using the STAR aligner.
############
configfile:"config.yaml"

gz_command="--readFilesCommand zcat" if config["gzipped"] else ""

path = config['path']
project = config['project']
organism = config['org']

# tsting for reading in a list pf files from a specific folder:
IDS, = glob_wildcards("rawData/{sample}.conc.R1.fastq.gz")

print(IDS)

rule all:
    input:
        expand('{project}/{organism}/star/bamFiles/{sample}.bam', sample = IDS),
        expand('{project}/{organism}/star/bamFiles/{sample}.bam.bai', sample = IDS),
#        expand("Analysis/{organism}/star/bwig/{sample}.bw", sample= IDS, organism = config['organism'])


rule map_star:
    input:
        R1='rawData/{IDS}.conc.R1.fastq.gz',
#
        index=expand("genome/{organisms}/starIndex/")
    output:
        bam='{project}/{organism}/star/bamFiles/{IDS}.bam',
        counts="{project}/{organism}/star/bamFiles/{IDS}.counts.tab",
        SJ = "{project}/{organism}/star/bamFiles/{IDS}SJ.out.tab",
    benchmark:
        "{project}/{organism}/benchmarks/{IDS}.run_rRNA_STAR.txt"
    params:
        prefix ="{project}/{organism}/star/bamFiles/{IDS}",
        gz_support=gz_command
    threads: 16
    shell:
        r'''
	mkdir -p {params.prefix}

####### for single-end reads samples data set
        STAR --runThreadN {threads} --genomeDir {input.index} --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {params.gz_support} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM {config[RAM]} --quantMode GeneCounts --outReadsUnmapped Fastx
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.prefix}ReadsPerGene.out.tab {output.counts}
        rmdir {params.prefix}
        '''

rule index:
    input:
        '{project}/{organism}/star/bamFiles/{IDS}.bam'
    output:
        '{project}/{organism}/star/bamFiles/{IDS}.bam.bai'
    shell:
        'samtools index {input}'

rule chrom_size:
    input:
        fastA = "genomes/{organism}.fa"
    output:
        fai = "genomes/{organism}.fa.fai",
        chromSize = "genome/{organism}.chromSize"
    shell:
        """
        samtools faidx {input.fastA}
        cut -f 1,2  {output.fai} > {output.chromSize}
        """

rule create_bigwig:
    input:
        bam = "{project}/{organism}/star/bamFiles/{IDS}.bam",
        bai = "{project}/{organism}/star/bamFiles/{IDS}.bam.bai"
    output:
        bw = "{project}/{organism}/star/bwig/{IDS}.bw"
    params:
        chromSize = "genomes/{organism}.chromSize",
        dir ="{project}/{organism}/star/bwig",
        prefix ="{project}/{organism}/star/bwig/{sample}"
    shell:
        """
        mkdir -p {params.dir}
        bam2wig.py  -i {input.bam} -s {params.chromSize} -o {params.prefix} &> {params.prefix}.log
        rm -f {wildcards.sample}.wig
        """

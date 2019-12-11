# the mpping file for either paired-end or songle-end mapping step using the STAR aligner.
############
configfile:"/fs/pool/pool-bcfngs/scripts/config.yaml"

gz_command="--readFilesCommand zcat" if config["gzipped"] else ""

path=config['path']
print("Directory of raw fastq reads:")
print(path)

# testing for reading in a list pf files from a specific folder:
IDS, = glob_wildcards("{sample}.conc.R1.fastq.gz")

print("Sample list:")
print(IDS)

project=config['project']
organism=config['org']
print("Project Number is:" ,project)
print("Mapping against genome", organism)

rule all:
    input:
        expand('{project}/{organism}/star/{sample}.bam', sample = IDS, organism = config['org'], project = config['project']),
        expand('{project}/{organism}/star/{sample}.bam.bai', sample = IDS, organism = config['org'], project = config['project']),
        expand("{project}/{organism}/bwig/{sample}.bw", sample= IDS, organism = config['org'], project = config['project'])

rule map_star:
    input:
        R1='{IDS}.conc.R1.fastq.gz'
    output:
        bam='{project}/{organism}/star/{IDS}.bam',
        counts="{project}/{organism}/star/{IDS}.counts.tab",
        SJ = "{project}/{organism}/star/{IDS}SJ.out.tab",
    benchmark:
        "{project}/{organism}/star/{IDS}.benchmarks.STAR.txt"
    params:
        prefix ="{project}/{organism}/star/{IDS}",
        gz_support=gz_command,
        index=expand("/fs/pool/pool-bcfngs/genomes/{organism}/starIndex/", organism = config['org'])
    threads: 16
    shell:
        r'''
	mkdir -p {params.prefix}
####### for paired-end reads samples data set
        STAR --runThreadN {threads} --genomeDir {params.index} --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {params.gz_support} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM {config[RAM]} --quantMode GeneCounts --outReadsUnmapped Fastx
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.prefix}ReadsPerGene.out.tab {output.counts}
        rmdir {params.prefix}
        '''


rule index:
    input:
        '{project}/{organism}/star/{IDS}.bam'
    output:
        '{project}/{organism}/star/{IDS}.bam.bai'
    shell:
        'samtools index {input}'

rule chrom_size:
    input:
        fastA = "/fs/pool/pool-bcfngs/genomes/{organism}.fa"
    output:
        fai = "/fs/pool/pool-bcfngs/genomes/{organism}.fa.fai",
        chromSize = "/fs/pool/pool-bcfngs/genomes/{organism}.chromSize"
    shell:
        """
        samtools faidx {input.fastA}
        cut -f 1,2  {output.fai} > {output.chromSize}
        """

rule create_bigwig:
    input:
        chromSize = "/fs/pool/pool-bcfngs/genomes/{organism}.chromSize",
        bam = "{project}/{organism}/star/{IDS}.bam",
        bai = "{project}/{organism}/star/{IDS}.bam.bai"
    output:
        bw = "{project}/{organism}/bwig/{IDS}.bw"
    params:
        dir ="{project}/{organism}/bwig",
        prefix ="{project}/{organism}/bwig/{IDS}",
        wig = "{project}/{organism}/bwig/{IDS}.wig"
    shell:
        """
        mkdir -p {params.dir}
        bam2wig.py  -i {input.bam} -s {input.chromSize} -o {params.prefix} &> {params.prefix}.log
        rm {params.wig}
        """

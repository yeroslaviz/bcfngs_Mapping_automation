# the mpping file for either paired-end or songle-end mapping step using the STAR aligner.
############
configfile:"/fs/pool/pool-bcfngs/scripts/config.yaml"

gz_command="--readFilesCommand zcat" if config["gzipped"] else ""

path=config['path']
print(path)

project=config['project']
organism=config['org']

# testing for reading in a list pf files from a specific folder:
IDS, = glob_wildcards("{sample}.conc.R1.fastq.gz")

print(IDS)

rule all:
    input:
        expand('/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{sample}.bam', sample = IDS, organism = config['org'], project = config['project']),
        expand('/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{sample}.bam.bai', sample = IDS, organism = config['org'], project = config['project']),
        expand("/fs/pool/pool-bcfngs/{project}/{organism}/star/bwig/{sample}.bw", sample= IDS, organism = config['org'], project = config['project'])


rule map_star:
    input:
        R1='{IDS}.conc.R1.fastq.gz',
        R2='{IDS}.conc.R2.fastq.gz',
    output:
        bam='/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}.bam',
        counts="/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}.counts.tab",
        SJ = "/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}SJ.out.tab",
    benchmark:
        "/fs/pool/pool-bcfngs/{project}/{organism}/benchmarks/{IDS}.run_rRNA_STAR.txt"
    params:
        prefix ="/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}",
        gz_support=gz_command,
        index=expand("/fs/pool/pool-bcfngs/genomes/{organism}/starIndex/", organism = config['org'])
    threads: 16
    shell:
        r'''
	mkdir -p {params.prefix}
####### for paired-end reads samples data set
        STAR --runThreadN {threads} --genomeDir {params.index} --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2} {params.gz_support} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM {config[RAM]} --quantMode GeneCounts --outReadsUnmapped Fastx
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.prefix}ReadsPerGene.out.tab {output.counts}
        rmdir {params.prefix}
        '''

rule index:
    input:
        '/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}.bam'
    output:
        '/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}.bam.bai'
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
        bam = "/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}.bam",
        bai = "/fs/pool/pool-bcfngs/{project}/{organism}/star/bamFiles/{IDS}.bam.bai"
    output:
        bw = "/fs/pool/pool-bcfngs/{project}/{organism}/star/bwig/{IDS}.bw"
    params:
        chromSize = "/fs/pool/pool-bcfngs/genomes/{organism}.chromSize",
        dir ="/fs/pool/pool-bcfngs/{project}/{organism}/star/bwig",
        prefix ="/fs/pool/pool-bcfngs/{project}/{organism}/star/bwig/{IDS}"
    shell:
        """
        mkdir -p {params.dir}
        bam2wig.py  -i {input.bam} -s {params.chromSize} -o {params.prefix} &> {params.prefix}.log
        rm -f *.wig
        """

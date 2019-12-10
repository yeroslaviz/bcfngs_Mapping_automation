
configfile: "config.yaml"

#print(config['organism'])
#print(config['organism']['Dmel'])
print(config['organism']['Dmel']['fasta'])
#print(config['organism']['Dmel'].keys())

rule all:
    input:
        expand("genome/{org}/starIndex/", org=config['organism']),
#        expand("genome/{org}/bwaIndex/", org=config['organism']),
#        expand("genome/{org}/bowtie2Index/", org=config['organism'])

#### get reference genomic data for the mapping (fasta and gtf files). Links are added in the config file
rule get_genome:
   output:
         fastA="genome/{org}.fa",
         gtf="genome/{org}.gtf"
   message:
          "getting {config[organism]} data - fastA and gtf files from Ensembl"
   params:
         fasta = lambda wildcards: config['organism'][wildcards.org]['fasta'],
         gtf   = lambda wildcards: config['organism'][wildcards.org]['gtf']
   shell:
         """
         wget -nc -O - {params.fasta} | gunzip -c - > {output.fastA}
         wget -nc -O - {params.gtf}   | gunzip -c - > {output.gtf}
         """

### Indexing the reference genome

rule star_index:
   input:
          fasta="genome/{org}.fa",
          gtf="genome/{org}.gtf"
   output:
          directory("genome/{org}/starIndex/")
   threads: 16
   params:
         prefix = lambda wildcards: "{org}".format(org=wildcards.org)
   shell:
          "mkdir -p {output} && "
          "STAR --runThreadN {threads} "
          "--outFileNamePrefix {output}{params.prefix} "
          "--runMode genomeGenerate "
          "--genomeDir {output} "
          "--limitGenomeGenerateRAM {config[RAM]} "
          "--genomeSAindexNbases {config[SAindex]} "
          "--genomeFastaFiles {input.fasta} "
          "--sjdbGTFfile {input.gtf} "
          "--sjdbOverhang 100"

rule bwa_index:
    input:
          fasta = "genome/{org}.fa"
    output:
          directory("genome/{org}/bwaIndex/")
    params:
        prefix = lambda wildcards: "{org}".format(org=wildcards.org)
    shell:
        "bwa index -b {config[blockSize]} -p {output}{params.prefix} {input.fasta} "


# rule bowtie2_index
rule bowtie2_index:
    input:
          fasta = "genome/{org}.fa"
    output:
          directory("genome/{org}/bowtie2Index/")
    params:
        prefix = lambda wildcards: "{org}".format(org=wildcards.org)
    message:
        "Indexing config['organism'] genome for bowtie2 Mapping"
#    log:
#        "log/bwa_index.log"
    shell:
        "bowtie2-build {input.fasta} {output}{params.prefix}"

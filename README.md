# Mapping automation

workflow for the bcfngs core facility at 

This is an internal how to for the workflow of the automation for the bcfngs.

## genome indexing

This steps needs to be run onyl once. **unless** new genomes are added. 

this is all run on the `hpcl5002`.

1. log-in to the server.
2. If not done before, activate the conda environment inside a `screen`, so that one can log out from the server and still run the workflow in the background. More information about the is in the github web site.
3. open the `config.yaml` file, e.g. ` vi /fs/pool/pool-bcfngs/scripts/config.yaml` and add the genome links from the ensembl FTP site in this format:

```
  Dps.3.0.45:
    fasta: "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/drosophila_pseudoobscura/dna/Drosophila_pseudoobscura.Dpse_3.0.dna.toplevel.fa.gz"
    gtf: "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gtf/drosophila_pseudoobscura/Drosophila_pseudoobscura.Dpse_3.0.45.gtf.gz"
```

4. run the command `snakemake -nps ../scripts/getGenome_IndexGenome.Snakefile` to check that the new genome is added to the workflow (the `-n` is a dry-run which just output the command without performiong it.)
5. run the same command without the `-n` - `snakemake -ps ../scripts/getGenome_IndexGenome.Snakefile`. This will download the fastA and gtf files from the given links and create the index files for them as well. 


## Mapping

Mapping woill be done with the `STAR` aligner for a specific genome either in a paired-end or singel-end format. 

### Changing parameters before mapping

The following settng must be adjusted before running the script:

```
######## Project Setting
# To create a separate project results folder for each run, change the name of the project here.
project: "P152"
# this path would point to the raw data in fastq format. This must be changed each time for the correct data to be analyzed.
path: "/fs/pool/pool-bcfngs/fastq_files/P146/P152_ChIP_Asx_5_6/conc.fastq"
# Which organism to map the data to
org: "Dme.BDGP6.22"
```

The `project` and the `path` should point to the current project which is being mapped.

the `org` must be the organism the data is mapped against. The term must be identical to the term in the list of organisms from the config files.

### running the `star` mapping step. 

first change the working directory to the path of the actual project on the pool-bcfngs. e.g.

```
cd /fs/pool/pool-bcfngs/fastq_files/P000/P000_Testrun/conc.fastq/
```

Next run the snakemake script first as a dry-run (`-n`) and if all seems well than for real. 

e.g. 

```
# for single-end samples
snakemake -nps /fs/pool/pool-bcfngs/scripts/Star.MappingQuant.SingleEndFastq.Snakefile -j 50
# for paired-end samples
snakemake -ps /fs/pool/pool-bcfngs/scripts/Star.MappingQuant.PairedEndFastq.Snakefile -j 50
```




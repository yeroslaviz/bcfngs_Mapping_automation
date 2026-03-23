import os
from pathlib import Path

configfile: "configs/genomes.yaml"

# Keep external downloads on the submit/login node. On some clusters
# compute nodes do not allow outbound internet access.
localrules: download_sources

GENOMES_ROOT = config["genomes_root"]
DEFAULT_THREADS = int(config.get("default_threads", 16))
DEFAULT_RAM = int(config.get("default_ram_bytes", 0))
INDEX_VERSIONS = config.get("index_versions", {})
STAR_VER = INDEX_VERSIONS.get("star", "unknown")
BWA_VER = INDEX_VERSIONS.get("bwa", "unknown")
BOWTIE2_VER = INDEX_VERSIONS.get("bowtie2", "unknown")
BWA_INDEX_DEFAULT_ALGO = str(config.get("bwa_index_algorithm", "is"))
BWA_INDEX_DEFAULT_MEM_MB = int(config.get("bwa_index_mem_mb", 16000))

ORGANISMS = config["organisms"]

ORG_KEYS = list(ORGANISMS.keys())
if "org_subset" in config and config["org_subset"]:
    requested = [o.strip() for o in str(config["org_subset"]).split(",") if o.strip()]
    missing = sorted(set(requested) - set(ORG_KEYS))
    if missing:
        raise ValueError(f"Unknown organism(s) in org_subset: {', '.join(missing)}")
    ORG_KEYS = requested


def _release_for_org(org_key: str) -> str:
    org = ORGANISMS[org_key]
    if org.get("source") == "ensemblgenomes":
        division = org.get("division")
        return str(config["releases"]["ensemblgenomes"][division])
    return str(config["releases"]["ensembl"])


def _source_root(org_key: str) -> str:
    org = ORGANISMS[org_key]
    if org.get("source") == "ensemblgenomes":
        division = org["division"]
        rel = _release_for_org(org_key)
        return f"https://ftp.ensemblgenomes.ebi.ac.uk/pub/{division}/release-{rel}"
    rel = _release_for_org(org_key)
    return f"https://ftp.ensembl.org/pub/release-{rel}"


def _normalize_ensemblgenomes_url(url: str) -> str:
    return url.replace(
        "https://ftp.ensemblgenomes.org/",
        "https://ftp.ensemblgenomes.ebi.ac.uk/",
        1,
    )


def _species_prefix(species: str) -> str:
    parts = species.split("_")
    if not parts:
        return species
    return "_".join([parts[0].capitalize(), *parts[1:]])


def fasta_url(wc):
    org = ORGANISMS[wc.org]
    if "fasta_url" in org:
        return _normalize_ensemblgenomes_url(org["fasta_url"])
    root = _source_root(wc.org)
    species = org["species"]
    assembly = org["assembly"]
    kind = org.get("fasta_kind", "dna.toplevel")
    # Example: .../fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    prefix = _species_prefix(species)
    return f"{root}/fasta/{species}/dna/{prefix}.{assembly}.{kind}.fa.gz"


def gtf_url(wc):
    org = ORGANISMS[wc.org]
    if "gtf_url" in org:
        return _normalize_ensemblgenomes_url(org["gtf_url"])
    root = _source_root(wc.org)
    species = org["species"]
    assembly = org["assembly"]
    rel = _release_for_org(wc.org)
    gtf_kind = org.get("gtf_kind", "gtf")
    prefix = _species_prefix(species)
    # Example: .../gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
    return f"{root}/gtf/{species}/{prefix}.{assembly}.{rel}.{gtf_kind}.gz"


def sources_dir(wc):
    org = ORGANISMS[wc]
    rel = _release_for_org(wc)
    src = org.get("source", "ensembl")
    if src == "ensemblgenomes":
        src = f"ensemblgenomes-{org['division']}"
    return os.path.join(GENOMES_ROOT, "sources", wc, f"{src}-release-{rel}")


def source_basename(org_key: str):
    assembly = ORGANISMS[org_key]["assembly"]
    return f"{org_key}.{assembly}"


def source_fasta_path(org_key: str):
    return os.path.join(sources_dir(org_key), f"{source_basename(org_key)}.fa")


def source_gtf_path(org_key: str):
    return os.path.join(sources_dir(org_key), f"{source_basename(org_key)}.gtf")


def source_fasta_path_for_wc(wc):
    assembly = ORGANISMS[wc.org]["assembly"]
    if wc.assembly != assembly:
        raise ValueError(f"Assembly mismatch for {wc.org}: expected {assembly}, got {wc.assembly}")
    return source_fasta_path(wc.org)


def source_gtf_path_for_wc(wc):
    assembly = ORGANISMS[wc.org]["assembly"]
    if wc.assembly != assembly:
        raise ValueError(f"Assembly mismatch for {wc.org}: expected {assembly}, got {wc.assembly}")
    return source_gtf_path(wc.org)


def tool_index_dirname(tool: str):
    return {
        "star": "starIndex",
        "bwa": "bwaIndex",
        "bowtie2": "bowtie2Index",
    }[tool]


def versioned_index_dir(tool: str, wc):
    ver = INDEX_VERSIONS.get(tool, "unknown")
    assembly = ORGANISMS[wc]["assembly"]
    return os.path.join(GENOMES_ROOT, wc, f"{assembly}-{tool_index_dirname(tool)}-{ver}")


def bwa_index_algorithm(wc):
    return str(ORGANISMS[wc.org].get("bwa_index_algorithm", BWA_INDEX_DEFAULT_ALGO))


def bwa_index_mem_mb(wc):
    return int(ORGANISMS[wc.org].get("bwa_index_mem_mb", BWA_INDEX_DEFAULT_MEM_MB))


rule all:
    input:
        expand(os.path.join(GENOMES_ROOT, "{org}.chromSize"), org=ORG_KEYS),
        [os.path.join(versioned_index_dir("star", org), "SA") for org in ORG_KEYS],
        [os.path.join(versioned_index_dir("bwa", org), f"{org}.bwt") for org in ORG_KEYS],
        [os.path.join(versioned_index_dir("bowtie2", org), f"{org}.1.bt2") for org in ORG_KEYS],


rule download_sources:
    output:
        fasta=os.path.join(GENOMES_ROOT, "sources", "{org}", "{src_rel}", "{org}.{assembly}.fa"),
        gtf=os.path.join(GENOMES_ROOT, "sources", "{org}", "{src_rel}", "{org}.{assembly}.gtf"),
    log:
        os.path.join(GENOMES_ROOT, "logs", "{org}", "{src_rel}", "download_sources.log"),
    conda:
        "../env/environment.yml"
    params:
        dir=lambda wc: sources_dir(wc.org),
        fasta_url=fasta_url,
        gtf_url=gtf_url,
    message:
        "Downloading FASTA/GTF for {wildcards.org}"
    shell:
        r"""
        mkdir -p "{params.dir}" "$(dirname "{log}")"
        : > "{log}"
        {{
          echo "Downloading FASTA: {params.fasta_url}"
          curl --retry 3 --retry-delay 5 --retry-connrefused -fsSL "{params.fasta_url}" | gunzip -c > "{output.fasta}"
          echo "Downloading GTF: {params.gtf_url}"
          curl --retry 3 --retry-delay 5 --retry-connrefused -fsSL "{params.gtf_url}"   | gunzip -c > "{output.gtf}"
        }} >> "{log}" 2>&1
        """


rule star_index:
    input:
        fasta=source_fasta_path_for_wc,
        gtf=source_gtf_path_for_wc,
    output:
        sa=os.path.join(GENOMES_ROOT, "{org}", f"{{assembly}}-starIndex-{STAR_VER}", "SA"),
    log:
        os.path.join(GENOMES_ROOT, "logs", "{org}", "star_index.log"),
    conda:
        "../env/environment.yml"
    threads: DEFAULT_THREADS
    params:
        outdir=lambda wc: os.path.join(GENOMES_ROOT, wc.org, f"{wc.assembly}-starIndex-{STAR_VER}"),
        sa_index=lambda wc: ORGANISMS[wc.org]["sa_index"],
        ram=lambda wc: DEFAULT_RAM,
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        : > "{log}"
        {{
          STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir "{params.outdir}" \
            --limitGenomeGenerateRAM {params.ram} \
            --genomeSAindexNbases {params.sa_index} \
            --genomeFastaFiles "{input.fasta}" \
            --sjdbGTFfile "{input.gtf}" \
            --sjdbOverhang 100
        }} >> "{log}" 2>&1
        """


rule bwa_index:
    input:
        fasta=source_fasta_path_for_wc,
    output:
        bwt=os.path.join(GENOMES_ROOT, "{org}", f"{{assembly}}-bwaIndex-{BWA_VER}", "{org}.bwt"),
    log:
        os.path.join(GENOMES_ROOT, "logs", "{org}", "bwa_index.log"),
    conda:
        "../env/environment.yml"
    threads: 1
    resources:
        mem_mb=bwa_index_mem_mb
    params:
        outdir=lambda wc: os.path.join(GENOMES_ROOT, wc.org, f"{wc.assembly}-bwaIndex-{BWA_VER}"),
        prefix=lambda wc: os.path.join(GENOMES_ROOT, wc.org, f"{wc.assembly}-bwaIndex-{BWA_VER}", wc.org),
        algo=bwa_index_algorithm,
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        : > "{log}"
        {{
          bwa index -a "{params.algo}" -p "{params.prefix}" "{input.fasta}"
        }} >> "{log}" 2>&1
        """


rule bowtie2_index:
    input:
        fasta=source_fasta_path_for_wc,
    output:
        bt2=os.path.join(GENOMES_ROOT, "{org}", f"{{assembly}}-bowtie2Index-{BOWTIE2_VER}", "{org}.1.bt2"),
    log:
        os.path.join(GENOMES_ROOT, "logs", "{org}", "bowtie2_index.log"),
    conda:
        "../env/environment.yml"
    threads: DEFAULT_THREADS
    params:
        outdir=lambda wc: os.path.join(GENOMES_ROOT, wc.org, f"{wc.assembly}-bowtie2Index-{BOWTIE2_VER}"),
        prefix=lambda wc: os.path.join(GENOMES_ROOT, wc.org, f"{wc.assembly}-bowtie2Index-{BOWTIE2_VER}", wc.org),
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        : > "{log}"
        {{
          bowtie2-build --threads {threads} "{input.fasta}" "{params.prefix}"
        }} >> "{log}" 2>&1
        """


rule chrom_size:
    input:
        fasta=lambda wc: source_fasta_path(wc.org),
    output:
        chrom=os.path.join(GENOMES_ROOT, "{org}.chromSize"),
    log:
        os.path.join(GENOMES_ROOT, "logs", "{org}", "chrom_size.log"),
    conda:
        "../env/environment.yml"
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        : > "{log}"
        {{
          samtools faidx "{input.fasta}"
          cut -f 1,2 "{input.fasta}.fai" > "{output.chrom}"
        }} >> "{log}" 2>&1
        """

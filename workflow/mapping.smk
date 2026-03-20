import os
from pathlib import Path


ALLOWED_MAPPERS = {"star", "bwa", "bowtie2"}

PROJECT = str(config.get("project", "P000")).strip()
if not PROJECT:
    raise ValueError("Missing config key: project")

ORGANISM = str(config.get("organism", "")).strip()
if not ORGANISM:
    raise ValueError("Missing config key: organism")
if "," in ORGANISM:
    raise ValueError("Only one organism is supported per mapping run.")

FASTQ_DIR = Path(str(config.get("fastq_dir", config.get("path", ".")))).expanduser()
if not FASTQ_DIR.exists() or not FASTQ_DIR.is_dir():
    raise ValueError(f"FASTQ directory not found or not a directory: {FASTQ_DIR}")

OUTDIR_RAW = str(config.get("outdir", ".")).strip()
if not OUTDIR_RAW:
    OUTDIR_RAW = "."
OUTDIR = Path(OUTDIR_RAW).expanduser()

GENOMES_ROOT = Path(str(config.get("genomes_root", "/fs/pool/pool-bcfngs/genomes"))).expanduser()
READ_MODE = str(config.get("read_mode", "paired")).strip().lower()
if READ_MODE not in {"paired", "single"}:
    raise ValueError(f"Invalid read_mode={READ_MODE!r}; expected paired or single.")

THREADS = int(config.get("threads", 32))
if THREADS < 1:
    raise ValueError("threads must be >= 1")

RAM_BYTES = int(config.get("ram_bytes", 168632718037))
if RAM_BYTES < 1:
    raise ValueError("ram_bytes must be >= 1")

MAPPERS_RAW = str(config.get("mappers", "star"))
MAPPERS = []
for item in [x.strip().lower() for x in MAPPERS_RAW.split(",") if x.strip()]:
    if item not in ALLOWED_MAPPERS:
        raise ValueError(f"Unsupported mapper: {item}")
    if item not in MAPPERS:
        MAPPERS.append(item)
if not MAPPERS:
    raise ValueError("No mapper selected.")

BW_BIN_SIZE = int(config.get("bw_bin_size", 10))
if BW_BIN_SIZE < 1:
    raise ValueError("bw_bin_size must be >= 1")
BW_NORMALIZATION = str(config.get("bw_normalization", "CPM")).strip()
if not BW_NORMALIZATION:
    raise ValueError("bw_normalization must not be empty")

FASTQ_PATTERN_DEFS = [
    {
        "name": "aviti",
        "glob": "*_R1.fastq.gz",
        "r1_suffix": "_R1.fastq.gz",
        "r2_suffix": "_R2.fastq.gz",
        "read_cmd": "zcat",
    },
    {
        "name": "aviti",
        "glob": "*_R1.fastq",
        "r1_suffix": "_R1.fastq",
        "r2_suffix": "_R2.fastq",
        "read_cmd": "cat",
    },
    {
        "name": "nova_seq",
        "glob": "*_R1_001.fastq.gz",
        "r1_suffix": "_R1_001.fastq.gz",
        "r2_suffix": "_R2_001.fastq.gz",
        "read_cmd": "zcat",
    },
    {
        "name": "nova_seq",
        "glob": "*_R1_001.fastq",
        "r1_suffix": "_R1_001.fastq",
        "r2_suffix": "_R2_001.fastq",
        "read_cmd": "cat",
    },
    {
        "name": "l1_pattern",
        "glob": "*_L1_R1.fastq.gz",
        "r1_suffix": "_L1_R1.fastq.gz",
        "r2_suffix": "_L1_R2.fastq.gz",
        "read_cmd": "zcat",
    },
    {
        "name": "l1_pattern",
        "glob": "*_L1_R1.fastq",
        "r1_suffix": "_L1_R1.fastq",
        "r2_suffix": "_L1_R2.fastq",
        "read_cmd": "cat",
    },
    {
        "name": "conc",
        "glob": "*.conc.R1.fastq.gz",
        "r1_suffix": ".conc.R1.fastq.gz",
        "r2_suffix": ".conc.R2.fastq.gz",
        "read_cmd": "zcat",
    },
    {
        "name": "conc",
        "glob": "*.conc.R1.fastq",
        "r1_suffix": ".conc.R1.fastq",
        "r2_suffix": ".conc.R2.fastq",
        "read_cmd": "cat",
    },
]


def detect_fastq_pattern():
    for pattern in FASTQ_PATTERN_DEFS:
        r1_files = sorted(FASTQ_DIR.glob(pattern["glob"]))
        if not r1_files:
            continue

        samples = []
        r1_by_sample = {}
        r2_by_sample = {}
        for r1_file in r1_files:
            name = r1_file.name
            if not name.endswith(pattern["r1_suffix"]):
                continue
            sample = name[: -len(pattern["r1_suffix"])]
            if not sample:
                continue
            samples.append(sample)
            r1_by_sample[sample] = str(r1_file)
            r2_by_sample[sample] = str(FASTQ_DIR / f"{sample}{pattern['r2_suffix']}")

        if samples:
            return pattern, sorted(samples), r1_by_sample, r2_by_sample

    raise ValueError(
        f"No recognizable FASTQ pattern found in {FASTQ_DIR}. "
        "Supported patterns: *_R1.fastq(.gz), *_R1_001.fastq(.gz), *_L1_R1.fastq(.gz), *.conc.R1.fastq(.gz)"
    )


FASTQ_PATTERN, SAMPLES, R1_FASTQ, R2_FASTQ = detect_fastq_pattern()
print(f"Detected FASTQ pattern: {FASTQ_PATTERN['name']} ({len(SAMPLES)} sample(s)) from {FASTQ_DIR}")

if READ_MODE == "paired":
    missing_r2 = [sample for sample in SAMPLES if not Path(R2_FASTQ[sample]).exists()]
    if missing_r2:
        preview = ", ".join(missing_r2[:5])
        extra = "" if len(missing_r2) <= 5 else f" (+{len(missing_r2) - 5} more)"
        raise ValueError(
            f"Paired-end mode selected, but missing R2 FASTQ files for samples: {preview}{extra}"
        )

STAR_INDEX = Path(str(config.get("star_index_dir", GENOMES_ROOT / ORGANISM / "starIndex"))).expanduser()
BWA_PREFIX = Path(str(config.get("bwa_index_prefix", GENOMES_ROOT / ORGANISM / "bwaIndex" / ORGANISM))).expanduser()
BOWTIE2_PREFIX = Path(str(config.get("bowtie2_index_prefix", GENOMES_ROOT / ORGANISM / "bowtie2Index" / ORGANISM))).expanduser()
CHROM_SIZE = GENOMES_ROOT / f"{ORGANISM}.chromSize"

missing_refs = []
if "star" in MAPPERS and not (STAR_INDEX / "SA").exists():
    missing_refs.append(str(STAR_INDEX / "SA"))
if "bwa" in MAPPERS and not Path(f"{BWA_PREFIX}.bwt").exists():
    missing_refs.append(f"{BWA_PREFIX}.bwt")
if "bowtie2" in MAPPERS and not Path(f"{BOWTIE2_PREFIX}.1.bt2").exists():
    missing_refs.append(f"{BOWTIE2_PREFIX}.1.bt2")
if not CHROM_SIZE.exists():
    missing_refs.append(str(CHROM_SIZE))
if missing_refs:
    raise ValueError("Missing required genome/index file(s):\n  - " + "\n  - ".join(missing_refs))

RUN_ROOT = OUTDIR / PROJECT / ORGANISM
RUN_ROOT_S = str(RUN_ROOT)
STAR_DIR = str(RUN_ROOT / "star")
BWA_DIR = str(RUN_ROOT / "bwa")
BOWTIE2_DIR = str(RUN_ROOT / "bowtie2")
LOG_DIR = str(RUN_ROOT / "logs")

ALL_OUTPUTS = (
    expand(os.path.join(RUN_ROOT_S, "{mapper}", "{sample}.bam"), mapper=MAPPERS, sample=SAMPLES)
    + expand(os.path.join(RUN_ROOT_S, "{mapper}", "{sample}.bam.bai"), mapper=MAPPERS, sample=SAMPLES)
    + expand(os.path.join(RUN_ROOT_S, "{mapper}", "bigwig", "{sample}.bw"), mapper=MAPPERS, sample=SAMPLES)
)
if "star" in MAPPERS:
    ALL_OUTPUTS += expand(os.path.join(STAR_DIR, "{sample}.counts.tab"), sample=SAMPLES)
    ALL_OUTPUTS += expand(os.path.join(STAR_DIR, "{sample}.SJ.out.tab"), sample=SAMPLES)


def r1_for_sample(wc):
    return R1_FASTQ[wc.sample]


def r2_for_sample(wc):
    return R2_FASTQ[wc.sample]


def star_input_reads(wc):
    if READ_MODE == "paired":
        return [R1_FASTQ[wc.sample], R2_FASTQ[wc.sample]]
    return [R1_FASTQ[wc.sample]]


rule all:
    input:
        ALL_OUTPUTS


rule map_star:
    input:
        reads=star_input_reads
    output:
        bam=os.path.join(STAR_DIR, "{sample}.bam"),
        counts=os.path.join(STAR_DIR, "{sample}.counts.tab"),
        sj=os.path.join(STAR_DIR, "{sample}.SJ.out.tab"),
    log:
        os.path.join(LOG_DIR, "star", "{sample}.log")
    threads: THREADS
    params:
        outdir=STAR_DIR,
        prefix=os.path.join(STAR_DIR, "{sample}.star."),
        index=str(STAR_INDEX),
        ram=RAM_BYTES,
        read_cmd=FASTQ_PATTERN["read_cmd"],
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        : > "{log}"
        {{
          rm -f "{params.prefix}"* "{output.bam}" "{output.counts}" "{output.sj}"

          STAR --runThreadN {threads} \
            --genomeDir "{params.index}" \
            --outFileNamePrefix "{params.prefix}" \
            --readFilesIn {input.reads} \
            --readFilesCommand {params.read_cmd} \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM {params.ram} \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx

          ln "{params.prefix}Aligned.sortedByCoord.out.bam" "{output.bam}"
          ln "{params.prefix}ReadsPerGene.out.tab" "{output.counts}"
          ln "{params.prefix}SJ.out.tab" "{output.sj}"
        }} >> "{log}" 2>&1
        """


rule map_bwa:
    input:
        r1=r1_for_sample,
    output:
        bam=os.path.join(BWA_DIR, "{sample}.bam"),
    log:
        os.path.join(LOG_DIR, "bwa", "{sample}.log")
    threads: THREADS
    params:
        outdir=BWA_DIR,
        index=str(BWA_PREFIX),
        paired=(READ_MODE == "paired"),
        r2=r2_for_sample,
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        : > "{log}"
        if [ "{params.paired}" = "True" ]; then
          bwa mem -t {threads} "{params.index}" "{input.r1}" "{params.r2}" 2>> "{log}" | \
            samtools sort -@ {threads} -o "{output.bam}" - >> "{log}" 2>&1
        else
          bwa mem -t {threads} "{params.index}" "{input.r1}" 2>> "{log}" | \
            samtools sort -@ {threads} -o "{output.bam}" - >> "{log}" 2>&1
        fi
        """


rule map_bowtie2:
    input:
        r1=r1_for_sample,
    output:
        bam=os.path.join(BOWTIE2_DIR, "{sample}.bam"),
    log:
        os.path.join(LOG_DIR, "bowtie2", "{sample}.log")
    threads: THREADS
    params:
        outdir=BOWTIE2_DIR,
        index=str(BOWTIE2_PREFIX),
        paired=(READ_MODE == "paired"),
        r2=r2_for_sample,
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname "{log}")"
        : > "{log}"
        if [ "{params.paired}" = "True" ]; then
          bowtie2 -p {threads} -x "{params.index}" -1 "{input.r1}" -2 "{params.r2}" 2>> "{log}" | \
            samtools sort -@ {threads} -o "{output.bam}" - >> "{log}" 2>&1
        else
          bowtie2 -p {threads} -x "{params.index}" -U "{input.r1}" 2>> "{log}" | \
            samtools sort -@ {threads} -o "{output.bam}" - >> "{log}" 2>&1
        fi
        """


rule index_bam:
    input:
        bam=os.path.join(RUN_ROOT_S, "{mapper}", "{sample}.bam")
    output:
        bai=os.path.join(RUN_ROOT_S, "{mapper}", "{sample}.bam.bai")
    threads: 1
    shell:
        r"""
        samtools index "{input.bam}"
        """


rule create_bigwig:
    input:
        chrom_size=str(CHROM_SIZE),
        bam=os.path.join(RUN_ROOT_S, "{mapper}", "{sample}.bam"),
        bai=os.path.join(RUN_ROOT_S, "{mapper}", "{sample}.bam.bai"),
    output:
        bw=os.path.join(RUN_ROOT_S, "{mapper}", "bigwig", "{sample}.bw")
    log:
        os.path.join(LOG_DIR, "bamCoverage", "{mapper}", "{sample}.log")
    threads: min(THREADS, 8)
    params:
        bw_dir=os.path.join(RUN_ROOT_S, "{mapper}", "bigwig"),
        bin_size=BW_BIN_SIZE,
        normalization=BW_NORMALIZATION,
    shell:
        r"""
        mkdir -p "{params.bw_dir}" "$(dirname "{log}")"
        : > "{log}"
        bamCoverage \
          --bam "{input.bam}" \
          --outFileName "{output.bw}" \
          --outFileFormat bigwig \
          --binSize {params.bin_size} \
          --normalizeUsing {params.normalization} \
          --numberOfProcessors {threads} >> "{log}" 2>&1
        """

# Mapping automation

Internal workflow for genome index management and NGS mapping at the bcfngs core facility.

## 1) Environment

```bash
./scripts/bcfngs env create
./scripts/bcfngs env check
```

## 2) Genomes and indices

Dry-run checks:

```bash
./scripts/bcfngs genomes plan
```

Build:

```bash
./scripts/bcfngs genomes build --cores 20
```

Subset:

```bash
./scripts/bcfngs genomes build --cores 20 --organism Homo_sapiens,Mus_musculus
```

## 3) Mapping workflow (STAR/BWA/Bowtie2)

The active mapping workflow is `workflow/mapping.smk`, controlled via:

- `./scripts/bcfngs mapping plan`
- `./scripts/bcfngs mapping run`

The old STAR-only Snakefiles were archived under `archive/legacy/mapping/`.

### Required config (`configs/run.yaml`)

```yaml
project: "P000"
path: "/fs/pool/pool-bcfngs/fastq_files/P000/P000_Testrun/conc.fastq/"
outdir: "/fs/pool/pool-bcfngs/mapping"
organism: "Caenorhabditis_elegans"
threads: 16
ram_bytes: 168632718037
```

### Run examples

```bash
# dry-run with defaults (paired-end, mapper=star)
./scripts/bcfngs mapping plan

# execute with defaults
./scripts/bcfngs mapping run

# run all mappers
./scripts/bcfngs mapping run --mapper star,bwa,bowtie2

# single-end mode
./scripts/bcfngs mapping run --single-end

# keep intermediary files
./scripts/bcfngs mapping run --keep-intermediary-files
```

### One-run overrides (without editing config)

```bash
./scripts/bcfngs mapping plan \
  --project P152 \
  --organism Drosophila_melanogaster \
  --fastq-dir /fs/pool/pool-bcfngs/fastq_files/P152/conc.fastq \
  --outdir /fs/pool/pool-bcfngs/mapping
```

### SLURM mode

```bash
./scripts/bcfngs mapping plan --executor slurm --jobs 20
./scripts/bcfngs mapping run  --executor slurm --jobs 20
```

Defaults for SLURM mode:

- partition: `p.hpcl8`
- mail user: `yeroslaviz@biochem.mpg.de`
- mail type: `ALL`
- time: `04:00:00`

### FASTQ pattern detection

The mapping workflow auto-detects input naming in this order:

1. `*_R1.fastq.gz`
2. `*_R1_001.fastq.gz`
3. `*_L1_R1.fastq.gz`
4. `*.conc.R1.fastq.gz`

In paired-end mode, matching `R2` files are required.

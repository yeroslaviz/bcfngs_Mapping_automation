# Changelog

## 2026-03-04

- Added unified mapping workflow: `workflow/mapping.smk`.
- Added `./scripts/bcfngs mapping plan` and `./scripts/bcfngs mapping run`.
- Added mapper selection via `--mapper star,bwa,bowtie2` (default `star`).
- Added unified read-mode selection via `--paired-end` (default) / `--single-end`.
- Added FASTQ auto-detection with ordered pattern matching:
  - `*_R1.fastq.gz`
  - `*_R1_001.fastq.gz`
  - `*_L1_R1.fastq.gz`
  - `*.conc.R1.fastq.gz`
- Added single-organism enforcement for mapping runs.
- Replaced `bam2wig.py` conversion with `deeptools` `bamCoverage`.
- Added `--keep-intermediary-files` (default behavior removes intermediaries).
- Added mapping output layout:
  - `<outdir>/<project>/<organism>/<mapper>/<sample>.bam`
  - `<outdir>/<project>/<organism>/<mapper>/<sample>.bam.bai`
  - `<outdir>/<project>/<organism>/<mapper>/bigwig/<sample>.bw`
- Added optional SLURM execution for mapping (`--executor slurm`) with per-rule resource defaults.
- Added mapping documentation page: `mapping.qmd`.
- Updated navigation and examples to reference mapping wrapper commands.
- Archived legacy STAR-only Snakefiles to `archive/legacy/mapping/`.
- Added `outdir` to `configs/run.yaml`.
- Added `deeptools` dependency in `env/environment.yml`.

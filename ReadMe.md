# Exome Pipe

## Introduction

This is a nextflow-based pipeline built to analyze exome data. There are two main pathways which the pipe will run:

- Unrelated samples

- Trio samples *(although technically there is no limit, there could be numerous related samples)*

## Structure

Nextflow script containing pre-processing and instructions for running sarek and Exomiser on FastQ data. 

Pipeline Flow:

**Stage 1**:

1. **cleanBeds**: *Cleans up the BED file for the sample group*

1. **renameFastQs**: *Groups fastQ files by pairs and renames them*

2. **runPipe**: *Runs [sarek pipeline](https://nf-co.re/sarek)*


**Stage 2** (differs for multisample and triosample):

3. **produceHpoString**: *Produces a short string from hpo file using Proband as key*

4. **produceExomiserYAML**: *Produces a YAML analysis file specifically for Proband*

5. **produceExomiserBatch**: *Produce a batch.txt file which lists the analysis YAML files (only in multiSample pipe)*

5. **runExomiser**: *Runs [exomiser](https://exomiser.github.io/Exomiser/general/) using batch analysis file*

## Use of pipe.nf

Set `alias nextflow=/efs/sam/bin/nextflow`

Then copy and paste [pipe.nf](pipe.nf) into Haggis (`vi pipe.nf` - `i` to insert and then `:wq` to save and quit)

`-c`: Configuration file

`--bed`: BED file to use with samples (is always `SureSelect_v6.bed`)

`--fastq`: Directory containing `.fastq` files

`--hpo`: Directory containing hpo files for unrelated samples || Hpo file for proband in trio sample

`--ped`: `.ped` file containing pedigree. 

`--pipe`: multiSample || trioSample

Example Commands:
```bash
nextflow run pipe.nf -profile slurm \
-c /efs/sam/configScripts/slurm.config \
--bed /efs/sam/Macrogen_HN00115050/SureSelect_v6.bed \
--fastq /efs/sam/Macrogen_HN00115050/fastq \
--hpo /efs/sam/Macrogen_HN00115050/hpo \
--ped /efs/sam/Macrogen_HN00115050/ped \
--pipe multiSample
```
OR
```bash
nextflow run pipe.nf -profile slurm \
-c /efs/sam/configScripts/slurm.config \
--bed /efs/sam/Macrogen_HN00115050/SureSelect_v6.bed \
--fastq /efs/sam/trio_example/fastq \
--hpo /efs/sam/trio_example/hpo/141641.hpo \
--ped /efs/sam/trio_example/ped \
--pipe trioSample

```

Runtime will be roughly 6-7 hours.

## Input Assumptions

- Samples are named in the format `proband_1.fastq.gz`
- `.bed` files contain `chr` prefix and 2 lines of unneccessary headers. 

## slurm.config

Configuration parameters for slurm when running `pipe.nf` and `nf-core/sarek`.

## DEPRECATED: pipe_multiSample.nf/pipe_trioSample.nf

Each contains the test pipelines for analyzing trio or multi-sample data. Replaced by `--pipe` tag in [pipe.nf](pipe.nf).

## DEPRECATED: exomiser-template.yml

Template for exomiser analysis.

## Ideas for future development

- Add in function which automatically uploads data to database
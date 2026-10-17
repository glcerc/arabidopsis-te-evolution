
# Arabidopsis TE Evolution

### Comparative Transposable Element Annotation and Evolutionary Dynamics in *Arabidopsis thaliana*

## Overview

This repository contains a reproducible bioinformatics workflow for investigating transposable element (TE) diversity, genomic distribution, and evolutionary dynamics in *Arabidopsis thaliana*.

The project initially focuses on accession **Co-4** and will be extended to comparative analyses across multiple accessions.

The workflow integrates de novo TE annotation, LTR retrotransposon classification, sequence identity analysis, and genomic visualization.

## Research Objectives

This project aims to:

1. Characterize the transposable element composition of *A. thaliana* accessions.
2. Compare the genomic abundance and distribution of major TE superfamilies.
3. Classify Copia and Gypsy LTR retrotransposons into evolutionary clades.
4. Investigate sequence divergence and approximate TE insertion histories.
5. Identify candidate accession-specific TE expansions.
6. Explore potential relationships between TE insertions and gene annotations.

## Initial Dataset

| Property | Description |
|---|---|
| Species | Arabidopsis thaliana |
| Accession | Co-4 |
| Assembly method | Flye |
| Genome assembly length | 136,966,697 bp |
| Assembly sequences | 139 |
| TE annotation | EDTA |
| LTR classification | TEsorter / REXdb plant |

The Co-4 assembly was generated during the Genome and Transcriptome Assembly practical at the University of Bern.

## Analysis Workflow

### 1. TE Annotation

EDTA was used to identify and classify transposable elements.

Parameters included:

- `--step all`
- `--sensitive 1`
- `--anno 1`
- `--wholeelement 1`
- `--species others`

An Arabidopsis CDS reference was used to reduce the misclassification of coding sequences.

### 2. Full-Length LTR Classification

Intact LTR retrotransposons identified by EDTA were analyzed using TEsorter and the REXdb plant database.

The analysis assigns Copia and Gypsy elements to evolutionary clades based on conserved protein domains.

### 3. LTR Sequence Identity

The similarity between the two LTRs of each intact retrotransposon was analyzed to investigate relative insertion histories.

### 4. Genome-Wide TE Distribution

TE superfamily occupancy and genomic distribution are investigated using EDTA annotations and R-based visualization.

### 5. Comparative TE Dynamics

Planned analyses include sequence divergence landscapes, clade abundance comparisons, and accession-specific TE insertion candidates.

## Preliminary Results: Co-4

### Genome-Wide TE Composition

EDTA reported that 15.31% of the Co-4 assembly was occupied by interspersed repetitive elements.

| TE superfamily | Genomic occupancy |
|---|---:|
| Gypsy | 4.68% |
| Helitron | 3.21% |
| Mutator | 3.17% |
| Copia | 1.33% |
| CACTA | 0.84% |
| L1 | 0.76% |

### Full-Length LTR Retrotransposons

A total of 264 intact LTR candidates were analyzed.

| Classification | Count |
|---|---:|
| Copia | 131 |
| Gypsy | 113 |
| Unclassified | 20 |

The most abundant classified clades were:

- Ale: 65
- Athila: 53
- Tork: 21
- Retand: 21

Although Gypsy elements occupied a larger proportion of the genome, Copia showed a slightly higher number of classified intact LTR elements.

These preliminary findings require further comparative analysis before evolutionary explanations can be established.

## Repository Structure

```text
arabidopsis-te-evolution/
├── config/
│   └── example.env
├── scripts/
├── figures/
├── results/
├── research/
├── README.md
└── .gitignore
```

## Configuration

The workflow was developed for the University of Bern HPC environment.

Users should create their own configuration file:

```bash
cp config/example.env config/local.env
```

Edit `config/local.env` to specify:

- `COURSEDIR`: Shared software and reference resources.
- `WORKDIR`: Working directory for analysis results.
- `GENOME`: Absolute path to the genome assembly FASTA.

The example configuration must be adapted to the local system before running the workflow.

### Running the Analysis

Submit jobs from the repository root:

```bash
sbatch scripts/01-run_EDTA.sh
sbatch scripts/02-run_TEsorter_full_length_LTRs.sh
sbatch scripts/03-run_LTR_identity.sh
sbatch scripts/04-run_TEsorter_library.sh
```

These commands illustrate the workflow. Each downstream analysis requires its corresponding input files to have been generated successfully. The jobs should not all be submitted simultaneously.

The scripts depend on course-provided Apptainer containers and reference resources. They are not yet a fully portable installation outside the course HPC environment.

## Research Questions

### RQ1: LTR Retrotransposon Evolution
Do particular Copia and Gypsy clades show evidence of distinct insertion histories across Arabidopsis accessions?

### RQ2: Accession-Specific TE Expansions
Are specific TE families enriched in Co-4 compared with other accessions?

### RQ3: TE Genomic Distribution
How do the genomic distributions of Gypsy, Copia, Helitron, and TIR elements vary across accessions?

### RQ4: TE–Gene Relationships
Are particular TE families enriched near protein-coding genes, including genes involved in stress responses?

### RQ5: Assembly and Annotation Effects
How do assembly continuity and annotation methods affect estimated TE abundance and genomic occupancy?

These questions are exploratory hypotheses rather than established findings.

## Analysis Status and Limitations

The initial EDTA job was interrupted when the course reservation ended. A subsequent job reused completed intermediate results.

The second EDTA run generated the main annotation files, summary tables, and evaluation reports. However, Slurm reported an out-of-memory termination near the end of the workflow.

The outputs are therefore treated as preliminary until their integrity and reproducibility have been independently verified.

## References

- Lian et al. (2024). *A pan-genome of 69 Arabidopsis thaliana accessions reveals a conserved genome structure throughout the global species range.*
- Ou et al. (2019). *Benchmarking transposable element annotation methods for creation of a streamlined, comprehensive pipeline.* Genome Biology.
- Zhang et al. (2022). *TEsorter: an accurate and fast method to classify LTR-retrotransposons in plant genomes.* Horticulture Research.

## Acknowledgements

This project originated from genome assembly and annotation practicals at the University of Bern.

The course provided software environments, reference datasets, and analysis scripts. Course-provided materials remain attributed to their original authors.

## Project Status

**Ongoing — exploratory research and comparative genomic analysis.**


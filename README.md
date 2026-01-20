# NAFLD Shotgun Metagenomes

Code and analysis scripts for **Testerman & Li et al., 2021** - a shotgun metagenomics study of Non-Alcoholic Fatty Liver Disease (NAFLD) patient samples.

## Overview

This repository contains the complete bioinformatics pipeline and statistical analysis code used to process 36 shotgun metagenomes from NAFLD patients. The workflow includes quality control, taxonomic profiling, functional analysis, and statistical comparisons between disease groups.

## Repository Structure

```
.
├── scripts/                    # Processing pipeline scripts (numbered by order)
│   ├── 01_human_removal.sh         # Remove human reads with KneadData
│   ├── 02_metaphlan_taxonomic.sh   # Taxonomic profiling with MetaPhlAn3
│   ├── 03_metaphlan_viruses.sh     # Viral profiling with MetaPhlAn3
│   ├── 04_concatenate_reads.sh     # Merge paired-end reads for HUMAnN3
│   ├── 05_humann_functional.sh     # Functional profiling with HUMAnN3
│   └── 06_humann_normalize.sh      # Normalize HUMAnN3 output
├── analysis/                   # R analysis and results
│   ├── R_phyloseq_microviz_ancom.md    # Main analysis (diversity, composition, ANCOM)
│   ├── Maaslin2_NAFLD_Metagenomes.md   # Differential abundance with Maaslin2
│   └── github_files/                    # Auto-generated figures from R markdown
└── figures/                    # Publication figures
    ├── main/                       # Main manuscript figures
    │   ├── Figure_1.png               # Alpha/beta diversity and barplots
    │   └── Figure_2.png               # PCA with taxa arrows and iris plots
    └── supplementary/              # Supplemental figures
        └── Supplemental_Figure_*.png
```

## Pipeline Overview

### Stage 1: Quality Control
**Script:** `01_human_removal.sh`
**Tool:** KneadData with Trimmomatic and FastQC
**Purpose:** Remove human-associated DNA sequences from raw paired-end FASTQ files

### Stage 2: Taxonomic Profiling
**Scripts:** `02_metaphlan_taxonomic.sh`, `03_metaphlan_viruses.sh`
**Tool:** MetaPhlAn 3 with Bowtie2
**Purpose:** Generate taxonomic abundance profiles at multiple levels (phylum to species) including viral profiling

### Stage 3: Functional Profiling
**Scripts:** `04_concatenate_reads.sh`, `05_humann_functional.sh`
**Tool:** HUMAnN 3
**Purpose:** Identify gene families and metabolic pathway abundances

### Stage 4: Normalization
**Script:** `06_humann_normalize.sh`
**Tool:** humann_renorm_table
**Purpose:** Convert raw counts to relative abundance

### Stage 5: Statistical Analysis
**Files:** `analysis/R_phyloseq_microviz_ancom.md`, `analysis/Maaslin2_NAFLD_Metagenomes.md`
**Tools:** R (phyloseq, microViz, ANCOM-BC, Maaslin2, patchwork)
**Analyses performed:**
- Alpha diversity (Shannon index)
- Beta diversity (Bray-Curtis NMDS, PCA)
- PERMANOVA and beta dispersion
- Compositional barplots
- Differential abundance testing

## Software Requirements

### Bioinformatics Tools
- [KneadData](https://huttenhower.sph.harvard.edu/kneaddata/) - Quality control and host read removal
- [MetaPhlAn 3](https://huttenhower.sph.harvard.edu/metaphlan/) - Taxonomic profiling
- [HUMAnN 3](https://huttenhower.sph.harvard.edu/humann/) - Functional profiling
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) - Sequence alignment

### R Packages
- phyloseq
- microViz
- ANCOM-BC (ancombc)
- Maaslin2
- patchwork
- vegan
- tidyverse

## Usage

### Running the Processing Pipeline

Scripts are numbered in order of execution. Each script processes samples listed in a `file_list.txt` file (not included - create based on your sample names).

```bash
# 1. Quality control - remove human reads
bash scripts/01_human_removal.sh

# 2. Taxonomic profiling
bash scripts/02_metaphlan_taxonomic.sh

# 3. Viral profiling (uses bowtie2 output from step 2)
bash scripts/03_metaphlan_viruses.sh

# 4. Concatenate reads for HUMAnN3
bash scripts/04_concatenate_reads.sh

# 5. Functional profiling
bash scripts/05_humann_functional.sh

# 6. Normalize output
bash scripts/06_humann_normalize.sh
```

### Expected Input Files

Scripts expect the following directory structure for input data:
- `Raw_Fastq/` - Raw paired-end FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`)
- `file_list.txt` - List of sample names (one per line)

### Running the R Analysis

The R analysis markdown files can be executed in RStudio or rendered using:

```r
rmarkdown::render("analysis/R_phyloseq_microviz_ancom.md")
```

## Data Availability

Raw sequencing data is available from [add accession number/repository link].

## Citation

If you use this code, please cite:

> Testerman T, Li Z, et al. (2021). [Paper title]. [Journal]. DOI: [add DOI]

## License

[Add license information]

## Contact

Todd Testerman - [add contact information]

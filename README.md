# Construction of RNA reference materials for improving the quality of transcriptomic data

RNA reference materials and the corresponding reference datasets that act as the “ground truth” of the measurement values are indispensable tools for assessing the reliability of RNA-seq in detecting intrinsically small biological differences in clinical settings such as those between molecular subtypes of diseases. However, constructing RNA reference datasets is challenging because of the incomparability of conventional “absolute” expression profiles across different batches, methods, or platforms. We recently proposed a ratio-based method for constructing reference datasets. The ratio for a gene is defined as the division of the expression levels between two sample groups and has demonstrated much better agreement than “absolute” data across multiple transcriptomic technologies and batches, resulting in the successful generation of omics-wide reference datasets with satisfied uncertainty.

In this protocol, we provide a step-by-step description of the procedures of establishing RNA reference materials and reference datasets, covering the following three stages: (1) reference materials, including material preparation, homogeneity testing, and stability testing; (2) ratio-based reference datasets, including characterization, uncertainty estimation, and orthogonal validation; and (3) applications, including construction of performance metrics, performing proficiency test, and diagnosing and correcting batch effects.

The protocol has been used for establishing the Quartet RNA reference materials and reference datasets (chinese-quartet.org) which have been approved as the first suite of certified RNA reference materials by China’s State Administration for Market Regulation as the First Class of National Reference Materials. The protocol can be utilized to establish and/or apply reference materials to improve RNA-seq data quality in diverse clinical utilities.

## About

Welcome to our repository dedicated to the construction of the RNA reference dataset. This repository contains the scripts utilized for this purpose, complemented by detailed instructions outlined in the accompanying article. Below, we provide comprehensive guidance on utilizing these scripts and interpreting the results effectively.

## Table of contents

- [Dependencies](#dependencies)
- [Setup](#setup)
- [Usage](#usage)
  * [RNA-seq](#RNA-seq)
  * [Examples](#examples)
    + [Starting with FASTQ files](#starting-with-fastq-files)
    + [Starting with BAM files](#starting-with-bam-files)
    + [Running prep and post separately](#running-prep-and-post-separately)
    + [Using the paired stats model](#using-the-paired-stats-model)
    + [Running the statistical model separately](#running-the-statistical-model-separately)
  * [Tips](#tips)
  * [All arguments](#all-arguments)
- [Output](#output)

## Dependencies

Tested on Ubuntu (20.04 LTS)
- Common
  - Python (3.10.14 and 2.7.18)
  - conda 24.3.0
    - R version 4.0.2 for standard R analysis
    - R version 4.3.3 for PVCA
    - Python (3.9.15 and 2.7.18)

- RNA-seq
  - Ballgown, version 2.20.0 (https://bioconductor.org/packages/release/bioc/html/ballgown.html)
  - Fastp, version 0.19.6 (https://github.com/OpenGene/fastp)
  - FastQ Screen, version 0.15.3 (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
  - FastQC, version 0.12.1 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - HISAT, version 2.1.0 (https://ccb.jhu.edu/software/hisat/index.shtml)
  - MultiQC version 1.8 (https://multiqc.info/)
  - Qualimap, version v2.3 (http://qualimap.conesalab.org/)
  - R, version 3.6.5 or higher (https://www.r-project.org/) and Rstudio (https://posit.co/downloads/) or any other integrated development environment
  - SAMtools, version 1.3.1 (https://github.com/samtools/samtools)
  - StringTie, version 1.3.4 (https://ccb.jhu.edu/software/stringtie/)


## Setup

If the required dependencies are already installed, then create a conda environment for Python and R dependencies:
```
# For RNA-seq upstream analysis
~/miniconda3/bin/conda create -n rna-seq-env python=3.7

# For PVCA only
~/miniconda3/bin/conda create -n pvca_env r bioconductor-pvca=1.42.0
```

And then run with:
```
# Activate the rna-seq-env conda environment
source ~/miniconda3/bin/activate ~/miniconda3/envs/rna-seq-env

# Install the requried tool packages
conda install -c bioconda fastp=0.19.6
conda install -c bioconda fastq-screen=0.12.0
conda install -c bioconda qualimap=2.3
conda install -c bioconda multiqc=1.8
conda install -c bioconda hisat2=2.1.0
conda install -c bioconda samtools=1.3.1
conda install -c bioconda stringtie=1.3.4
conda install -c bioconda bioconductor-ballgown=2.20.0

# Deactivate the rna-seq-env conda environment
conda deactivate
```

## Usage
### RNA-seq
▲CRITICAL Detailed procedures are provided in the Supplementary methods. All the command examples use the raw fastq data of sample D5_1 from batch R_BGI_L3_B1 in this protocol.
#### Alignment and gene quantification 
1.	Use fastp to remove adapter sequences from fastq files. For a detailed description of fastp, visit https://github.com/mdshw5/fastqp
2.	Perform RNA-seq alignment and gene-level quantification using HISAT, StringTie, and Ballgown tools. Specifically, use HISAT to map the trimmed reads to the reference genome, such as GRCh38. Then use StringTie to assemble and estimate gene abundances for each sample. Next, use Ballgown to create read-count and FPKM (Fragments Per Kilobase of transcript per Million mapped reads) matrices. Pertea et al. provided a detailed protocol for these tools.
3.	Add a value of 0.01 to the FPKM value of each gene and perform log2-transformation.
4.	Use expression profiles based on detected genes for further analysis. In our study, a gene is considered detectable (expressed) in a biological group within a batch if ≥ 3 reads were mapped onto it in at least two of the three replicates.

▲CRITICAL Ballgown can process a single sample or merge the expression profiles of multiple samples.

#### Quality control analysis
1.	Perform quality control on sequence data in fastq files using FastQC. For a detailed description of FastQC, visit https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
2.	Use FastQC on the fastq files generated by fastp to verify the success of adapter removal.
3.	Use FastQ Screen on the fastq files generated by fastp to detect potential contamination with other species, junction primers, etc. For a detailed description of FastQ Screen, visit https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
4.	Use Qualimap2 to calculate the quality of the reads that are mapped to different genomic regions. For a detailed description of Qualimap2, visit http://qualimap.conesalab.org/.

▲CRITICAL Qualimap takes longer to run and requires more computational resources, whereas the RNA-seq QC module may not run successfully due to version issues.

5.	Use MultiQC to merge the QC results from FastQC, FastQ Screen and Qualimap. For a detailed description of MultiQC, visit https://multiqc.info/.

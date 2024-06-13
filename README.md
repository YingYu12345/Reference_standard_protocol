# Construction of RNA reference materials for improving the quality of transcriptomic data

RNA reference materials and the corresponding reference datasets that act as the “ground truth” of the measurement values are indispensable tools for assessing the reliability of RNA-seq in detecting intrinsically small biological differences in clinical settings, such as those between molecular subtypes of diseases. However, constructing RNA reference datasets is challenging because conventional “absolute” expression profiles are incomparable across different batches, methods, or platforms. We recently proposed a ratio-based method for constructing reference datasets. The ratio for a gene is defined as the division of the expression levels between two sample groups and has demonstrated much better agreement than “absolute” data across multiple transcriptomic technologies and batches, resulting in the successful generation of omics-wide reference datasets with satisfied uncertainty.

In this protocol, we provide a step-by-step description of the procedures for establishing RNA reference materials and reference datasets, covering the following three stages: (1) **reference materials**, including material preparation, homogeneity testing, and stability testing; (2) **ratio-based reference datasets**, including characterization, uncertainty estimation, and orthogonal validation; and (3) **applications**, including definition of performance metrics, performing proficiency test, and diagnosing and correcting batch effects. 

The protocol has been used to establish the Quartet RNA reference materials and reference datasets ([chinese-quartet.org](https://chinese-quartet.org)) that have been approved by China's State Administration for Market Regulation as the first suite of certified RNA reference materials as the First Class of National Reference Materials. The protocol can be utilized to establish and/or apply reference materials to improve RNA-seq data quality in diverse clinical utilities.


## About
Chinese Quartet Multi-Omics Project Official Website: [chinese-quartet.org](https://chinese-quartet.org)

Welcome to our repository dedicated to the construction of the RNA reference dataset. This repository contains the scripts for this purpose, complemented by detailed instructions outlined in the accompanying article. Below, we provide comprehensive guidance on utilizing these scripts and interpreting the results effectively.

![fig1_overview_of_the_protocol](https://github.com/YingYu12345/Reference_standard_protocol/assets/72855092/c1e24325-76a0-48e5-84d7-fbfd48deecdd)
  **Overview of the protocol**

Establishing RNA materials and reference datasets involves three stages: (a) reference materials, (b) ratio-based reference datasets, and (c) applications.

## Table of contents

- [Dependencies](#dependencies)
- [Setup](#setup)
- [Usage](#usage)
  * [RNA-seq](#rna-seq)
  * [Homogeneity testing](#homogeneity-testing)
  * [Stability testing](#stability-testing)
  * [Characterization](#characterization)
  * [Measurement uncertainty estimation](#measurement-uncertainty-estimation)
  * [Construction of performance metrics](#construction-of-performance-metrics)
  * [Diagnosis and correction of batch effects](#diagnosis-and-correction-of-batch-effects)
- [Output](#output)
  * [Homogeneity](#homogeneity)
  * [Stability](#stability)
  * [Charac](#charac)
  * [Uncertainty](#uncertainty)
  * [Performance metrics](#performance-metrics)
  * [Batch effects](#batch-effects)
- [Reference](#reference)
- [Contact](#contact)

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

If the required common dependencies are already installed, then create a conda environment for Python and R dependencies:
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

> The following provides detailed instructions only for using the involved R scripts. For step-by-step details, please refer to the accompanying protocol article.

### RNA-seq
▲CRITICAL Detailed procedures are provided in the Supplementary Methods section. All the command examples use the raw fastq data of sample D5_1 from batch R_BGI_L3_B1 in this protocol.
#### Alignment and gene quantification 
1.	Use fastp to remove adapter sequences from fastq files. For a detailed description of fastp, visit https://github.com/OpenGene/fastp.
2.	Perform RNA-seq alignment and gene-level quantification using HISAT, StringTie, and Ballgown tools. Specifically, use HISAT to map the trimmed reads to the reference genome, such as GRCh38. Then use StringTie to assemble and estimate gene abundances for each sample. Next, use Ballgown to create read-count and FPKM (Fragments Per Kilobase of transcript per Million mapped reads) matrices. [Pertea et al.](https://www.nature.com/articles/nprot.2016.095) provided a detailed protocol for these tools.
3.	Add a value of 0.01 to the FPKM value of each gene and perform log2-transformation.
4.	Use expression profiles based on detected genes for further analysis. In our study, a gene is considered detectable (expressed) in a biological group within a batch if ≥ 3 reads were mapped onto it in at least two of the three replicates.

▲CRITICAL Ballgown can process a single sample or merge the expression profiles of multiple samples.

#### Quality control analysis
1.	Perform quality control on sequence data in fastq files using FastQC. For a detailed description of FastQC, visit https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
2.	Use FastQC on the fastq files generated by fastp to verify the success of adapter removal.
3.	Use FastQ Screen on the fastq files generated by fastp to detect potential contamination with other species, junction primers, etc. For a detailed description of FastQ Screen, visit https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
4.	Use Qualimap2 to calculate the quality of the reads mapped to different genomic regions. For a detailed description of Qualimap2, visit http://qualimap.conesalab.org/.

▲CRITICAL Qualimap takes longer to run and requires more computational resources, whereas the RNA-seq QC module may not run successfully due to version issues.

5.	Use MultiQC to merge the QC results from FastQC, FastQ Screen, and Qualimap. For a detailed description of MultiQC, visit https://multiqc.info/.

### Homogeneity testing
1.	Perform statistical analysis (e.g., ANOVA) to compare the variation between the expression profiles of genes within-unit and those between-unit groups. The expression profiles are obtained from RNA-seq analysis.
2.	Calculate the adjusted p values with the "fdr" method (e.g., "fdr" for Benjamini-Hochberg). A gene is considered homogeneous if the FDR-adjust ANOVA-based p >0.05, indicating no significant difference between within-unit and between-unit groups.

▲CRITICAL The homogeneous analysis used here is relatively loose. Users can choose a more robust test of homogeneity.

3.	Calculate the proportion of genes meeting the homogeneity criterion. If a majority of genes (e.g.,>95%) show homogeneity, the reference material is recognized as homogeneous.

Run the script homogeneity_assess.R to execute:

```
Rscript path-to/homogeneity_assess.R -i path-to/example_expr_homo_log2.csv -m path-to/metadata_homo.csv -o path-to/
```
▲CRITICAL Two files will be obtained after successfully implementing the script, including: (1) the proportion of genes that meet the homogeneity criterion (homogeneity_Ftest_summary_(date).csv); (2) the results of homogeneity of each gene under estimation (homogeneity_Ftest_detail_(date).csv). 

**PS: `path-to/` should be replaced with the actual file path, which applies to all instances below.**

### Stability testing
Calculate the ratio of RIN values using one reference material as the denominator. Perform regression analysis for each gene using lm() function to estimate the observed slope and uncertainty of slope of RIN values across time points. If the observed slope is smaller than the uncertainty of the slope multiplied by the confidence level, the materials are stable. 

Run the script stability_assess.R:
```
Rscript path-to/stability_assess.R -i path-to/example_expr_stability_RIN.csv -m path-to/metadata_stability_RIN.csv -o path-to/
```
▲CRITICAL One file named stability_lm_detail_(date).csv)will be obtained after successfully implementing the script, illustrating the results of stability of each gene or RNA RIN values under estimation. 

### Characterization
![fig4_workflow_of_characterization_of_reference_datasets](https://github.com/YingYu12345/Reference_standard_protocol/assets/72855092/bff73c7f-a111-443a-b6b7-bd21fd680595)
#### Data quality control based on expression profiles
1. Perform principal component analysis (PCA) to visualize the major sources of variation of high-dimensional data in each batch. Use prcomp() function in R.

2.	Calculate signal-to-noise ratio (SNR). SNR is a measure often expressed in decibels that quantifies the strength of a signal relative to the background noise. For RNA-seq data, the “signal” is the average distance that represents inherent “differences” among various biological sample groups, whereas the “noise” is the average distance among technical replicates within the same sample group in a space with reduced dimensionality. Generally, a lower SNR value suggests less discriminative power and vice versa. In our research, we have established the SNR threshold of 12. 

Run the script SNR_qc.R:
```
Rscript path-to/SNR_qc.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/
```
▲CRITICAL Three files will be obtained after successfully implementing the script, including (1) the SNR value of the expression profiles of each batch (SNR_perBatch_(date).csv); (2) scatter plots of the PCA results with the SNR value for each batch (SNR_(batch)_(date).pdf); and (3) a scatter plot of the PCA results with SNR value across all batches (SNR_allBatch_(date).pdf).

#### Identification of high-confidence detected genes
1.	Define expression criteria for considering a gene as expressed in a sample. A gene is considered expressed in a library in each batch if more than three reads were mapped to it in at least two of the three replicates. 
2.	Calculate the number of genes detected in each reference RNA sample group. 
3.	Conduct the detected genes across all samples. If a gene is detected beyond two of three replicates in a sample group of reference material, it is considered expressed in that sample group. 

Run the script char_detect.R to execute:

```
Rscript path-to/char_detect.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -o path-to/
```
▲CRITICAL A file named detect_genelist_(date).csv will be generated upon successful script implementation. The file contains a list of genes classified as detected genes in each batch within each reference RNA sample group.

#### Characterization of reference datasets at the ratio level 
1.	Identify the detectable genes across the two groups of each sample pair.
2.	Determine the ratio-based expressions based on detectable genes as follows:
a)	Obtain ratio-based expression data for each batch, calculated per-gene basis.
b)	Use log2FPKM values to calculate these ratio-based expressions.
c)	For each gene, first compute the mean of the expression profiles of replicates of one group of the reference materials.
d)	Subtract this mean from the log2FPKM values of that gene in each batch to obtain the ratio-based expressions.
3.	To improve the reliability of the reference values, select genes that meet the threshold of p < 0.05 in each sample pair using the limma method. 

Run the script deg_limma.R using the following command: 
```
Rscript path-to/deg_limma.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata char.csv -o path-to/
```
▲CRITICAL A file named DEG_limma_(date).csv will be generated upon successful script implementation. The file includes log2-transformed fold change (log2FC), averaged expression across both sample groups, t values, p values, FDR-adjusted p values, and differentially expressed gene classification in each sample pair. 

4.	Perform the Shapiro-Wilk test to check the normality of the ratio-based expressions. Ensure that data comply with normal distribution for statistical tests. 

Run the script normality.R using the following command: 
```
Rscript path-to/normality.R -i path-to/DEG_limma_(date).csv -p path-to/batch_included.csv -o path-to/
```
▲CRITICAL A file named normality_genelist_(date).csv will be generated upon successful script implementation. The file contains p values of the Shapiro-Wilk test based on the ratio-based expression of each gene in the reference datasets. 

5.	Generate reference datasets at the ratio level for each pair of reference materials in the format of a mean by summarizing high-quality RNA-seq datasets. 

Run the script char_ratio.R using the following command: 
```
Rscript path-to/char_ratio.R -i path-to/example_expr_multibatch_count.csv -m path-to/metadata.csv -g path-to/detect_genelist_(date).csv -d path-to/DEG_limma_(date).csv -o path-to/
```
▲CRITICAL A file named ref_expr_(date).csv will be generated upon successful script implementation. The file contains ratio-based expressions of reference datasets, providing information on FC and the median of limma-based p-value for each gene.

#### Identification of reference DEGs
Identify reference DEGs across high-quality batches. A gene is considered as a reference DEG between two sample groups if it is consistently categorized as an up- or down-regulated gene in more than a certain percentage of high-quality batches. In our study, we employed an arbitrary cutoff of 1/3. 

Run the script:
```
Rscript path-to/char_refDEG.R -d path-to/DEG_limma_(date).csv -r ref_expr_(date).csv -o path-to/
```
▲CRITICAL A file named RefData_DEGs_(date).csv will be generated upon successful script implementation. This file includes a list of genes classified as reference DEGs for each sample pair. Additionally, it lists the average FC and median p values across high-quality batches for each gene.

### Measurement uncertainty estimation
▲CRITICAL Follow the GUM (Guide to the Expression of Uncertainty in Measurement) method.

1.	Characterization uncertainty. Characterization uncertainty refers to the uncertainty associated with the characterization of reference datasets. It is calculated based on high-quality RNA-seq datasets. Specifically, for each gene generated in the characterization at the ratio level, uncertainty is estimated using fold changes (log2-transformed) from high-quality datasets. 

Run the script: 
```
Rscript path-to/uchar.R -i path-to/DEG_limma.csv -p path-to/passbatch.csv -r path-to/ref_expr_(date).csv -o path-to/
```
▲CRITICAL Upon successful script implementation, one file named  Refdata_uchar_(date).csv will be generated, containing the characterization uncertainty of each gene in the reference datasets.

2.	Inhomogeneity uncertainty. Inhomogeneity uncertainty arises when the material being measured is not uniform throughout. It is calculated using expression profiles. 

Run the script:
```
Rscript path-to/ubb.R -i path-to/example_expr_homo_log2.csv -m path-to/metadata_homo.csv -r path-to/ref_expr_(date).csv -o path-to/
```

▲CRITICAL Upon successful script implementation, one file named  Refdata_ubb_(date).csv will be generated, containing the characterization uncertainty of each gene in the reference datasets.

3.	Instability uncertainty. Instability uncertainty accounts for variations over time due to environmental conditions or other factors. It can be calculated using RNA-seq datasets or RIN values to illustrate overall RNA integrity. Our study employed linear regression profiles of RIN to represent overall instability. 

Run the script:
```
Rscript path-to/us.R -i path-to/example_expr_stability_RIN.csv -m path-to/metadata_stability_RIN.csv -r path-to/ref_expr_(date).csv -t R -o path-to/
```
▲CRITICAL Upon successful script implementation, one file named Refdata_us_(date).csv will be generated, containing the instability uncertainty of each gene in the reference datasets.

4.	Combined uncertainty. Combined uncertainty accounts for all sources affecting a measurement, including characterization uncertainty, inhomogeneity uncertainty, and instability uncertainty. 

Run the script: 
```
Rscript path-to/uc.R -c Refdata_uchar_(date).csv -b Refdata_ubb_(date).csv -s Refdata_us_(date).csv -r path-to/ref_expr_(date).csv -o path-to/
```
▲CRITICAL A file named Refdata_uc_(date).csv will be generated upon successful script implementation. These files contain the combined uncertainty of each gene in the reference datasets.

5.	Extended uncertainty. Extended uncertainty provides a confidence interval around the measured value. Apply a coverage factor (typically k = 2 for a 95% confidence level). Multiply the combined uncertainty by k to obtain the extended uncertainty. 

Run the script: 
```
Rscript path-to/U.R -d Refdata_uc_(date).csv -r path-to/ref_expr_(date).csv -k 2 -o path-to/
```
▲CRITICAL A file named Refdata_U_(date).csv will be generated upon successful script implementation. These files contain the extended uncertainty of each gene in the reference datasets.

### Construction of performance metrics
1.	Calculate relative correlation with reference datasets (RC) metrics. RC is calculated based on the Pearson correlation coefficient between the ratio-based expression levels of a dataset for a given pair of sample groups and the corresponding reference fold-change values (transformed to log2-scale). It is called the “relative correlation with reference datasets” metric, representing the numerical consistency of the ratio-based expression profiles. To improve reliability, calculate the mean of the three replicates of each sample group before performing ratio-based expression analysis. 

Run the script:
```
Rscript path-to/qc_rc.R -i path-to/example_log2_limma_data.csv -r path-to/example_ref_data.csv -o path-to/
```
▲CRITICAL A file named RC_each_batch_(date).csv will be generated upon successful script implementation. This file contains the RC values between each batch's test and reference datasets.

2.	Calculate the Matthews Correlation Coefficient (MCC) of reference DEGs. MCC measures the consistency of DEGs detected from a dataset for a given pair of samples with those from the reference DEGs, or “MCC of DEGs”, as follows:
a)	Repeat Step 86 to identify DEGs and non-DEGs of a test dataset.
b)	Obtain reference DEGs and non-DEGs from Steps 87-88. They can be true positive and true negative sets. 
c)	Compare DEGs and non-DEGs with reference DEGs and non-DEGs. Compute the number of True Positive (TP), True Negative (TN), False Positive (FP) and False Negative (FN).
d)	Calculate MCC using TP, TN, FP, and FN values.

Run the script: 
```
Rscript path-to/qc_mcc.R -i path-to/example_log2_limma_data.csv -r path-to/example_ref_data.csv -o path-to/
```
▲CRITICAL A file named MCC_each_batch_(date).csv will be generated upon successful script implementation. This file contains the TP, TN, FP, FN, and MCC values between each batch's test and reference datasets.

3.	Calculate the root mean square error (RMSE) using fold-changes between a test dataset for a given pair of samples and the corresponding ratio-based reference datasets, representing the average distances of ratio-based expression profiles. Fold-changes were transformed using log2 scaling. Use rmse() function Metrics package for implementation. 
Run the script: 
```
Rscript path-to/qc_rmse.R -i path-to/example_log2_limma_data.csv -r path-to/example_ref_data.csv -o path-to/
```
▲CRITICAL A file named RMSE_each_batch_(date).csv will be generated upon successful script implementation. This file contains the RMSE value between test and reference datasets for each batch.

### Diagnosis and correction of batch effects
1.	Use the principal variance component analysis (PVCA) method to assess whether batch effects are present in the dataset. When combined with bar plots, PVCA can quantify and visualize the proportion of variations attributed to experimental effects, including batch effects. 

Run the script: 
```
source ~/miniconda3/bin/activate ~/miniconda3/envs/pvca_env

Rscript path-to/pvca.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/
```
! CAUTION PVCA installation requires more dependent packages and may conflict with other R package versions. Creating a conda environment dedicated to PVCA analysis was recommended.

▲CRITICAL Two files will be obtained after successfully implementing the script, including (1) PVCA value of the expression profiles of multiple batches (PVCA_(date).csv) and (2) a bar plot of the PVCA results of the expression profiles of multiple batches (PVCA_plot_(date).pdf).

2.	A reference-material-based ratio method should be applied for batch-effect correction if batch effects are detected or suspected. Specifically, ratio-based expression data are obtained within each batch on a gene-by-gene basis. For each gene, calculate the mean of expression profiles (e.g. log2FPKM) of replicates of reference sample(s), and then subtract it from the log2FPKM values of that gene in each study. 
Run the script: 
```
Rscript path-to/ratio.R -i path-to/example_expr_multibatch_log2.csv -m path-to/metadata.csv -o path-to/
```
▲CRITICAL Samples that will be used as the denominator for performing the ratio-based method shall be annotated in the metadata.csv file. 

▲CRITICAL This step is not obligatory and depends on the preceding diagnostic results. We recommend to use log2-transformed FPKM data in this step. A file named RNA_expression_ratio_(date).csv will be obtained after successfully implementing the script, which contains the ratio-based expression profiles of each gene. 

## Output

### Homogeneity
(1) **homogeneity_Ftest_summary_(date).csv**: This file contains the proportion of genes that meet the homogeneity criterion.

The CSV file has the following columns:

1. **sample**: The sample ID (e.g., "D5").
2. **length_sig**: The number of genes that meet the homogeneity criterion (i.e., significant genes with FDR < 0.05).
3. **length_total**: The total number of genes tested.
4. **length_ratio**: The proportion of significant genes (`length_sig` / `length_total`).

The data in this file can be used to quickly assess the homogeneity of gene expression values across different samples. The columns provide insight into the number of significant genes (those that meet the homogeneity criterion) and their proportion concerning the total number of genes tested. This summary can help identify samples with a higher proportion of homogeneous gene expression patterns.

For more detailed information on individual gene results, refer to the corresponding `homogeneity_Ftest_detail_(date).csv` file.

(2) **homogeneity_Ftest_detail_(date).csv**: This file provides the results of the homogeneity test for each gene under estimation.

The CSV file has the following columns:

1. **sample**: The sample ID (e.g., "D5").
2. **gene**: The gene identifier (e.g., "ENSG00000239951").
3. **anova_p**: The p-value from the ANOVA test.
4. **mean_between**: The mean expression level between different groups.
5. **mean_intra**: The mean expression level is the same group.
6. **fdr**: The false discovery rate adjusted for multiple testing. A lower FDR indicates a higher confidence in the significance of the test results.

The data in this file can be used to assess the homogeneity of gene expression values between and within groups in different samples. The ANOVA p-values help determine the significance of differences, while the mean values and FDR provide additional context for the results.

### Stability
**stability_lm_detail_(date).csv**: This file llustrates the results of stability of each gene or RNA RIN values under estimation. 

The CSV file has the following columns:

1. **compare**: The comparison between sample pairs (e.g., "D6/D5").
2. **meanFC**: The mean fold change in RIN values.
3. **b1_slope**: The slope of the regression line (b1).
4. **b0_intercept**: The intercept of the regression line (b0).
5. **n**: The number of data points used in the regression analysis.
6. **s**: The standard error of the regression.
7. **sb1**: The standard error of the slope (b1).
8. **t95_n**: The critical value of the t-distribution at 95% confidence for n-2 degrees of freedom.
9. **sig**: A boolean indicating whether the slope (b1) is significantly different from zero (TRUE or FALSE).

The data in this file can be used to assess the stability of gene expression levels (take RIN values as an example) over time for different sample pairs. The columns provide detailed information about the regression analysis, including the mean fold change, the slope and intercept of the regression line, and the significance of the slope.

### Charac
#### the output of SNR_qc.R:
(1) **SNR_perBatch_(date).csv**: This file contains the SNR value of the expression profiles of each batch.

1. **batch_id**: Identifier for each batch. It includes information about the sequencing platform, library, and batch number.
2. **snr**: The calculated Signal-to-Noise Ratio for the each batch.

This file provides a summary of the Signal-to-Noise Ratio (SNR) analysis for different batches.

(2) **SNR_perBatch_(date).pdf**: This file provides scatter plots of the PCA results with the SNR value for each batch.

(3) **SNR_allBatch_(date).pdf**: This file provides a scatter plot of the PCA results with SNR value across all batches.

#### the output of char_detect.R:

**detect_genelist_(date).csv**: The file contains a list of genes classified as detected genes in each batch within each reference RNA sample group.

1. **gene**: This column lists the Ensembl gene identifiers for the genes that were detected.
2. **freq**: This column indicates the frequency of detection for each gene across different batches. A frequency of 10 suggests that the gene was detected in all batches analyzed.
3. **sample**: This column specifies the sample identifier in which the gene was detected. 

#### the output of deg_limma.R:
**DEG_limma_(date).csv**: The file includes log2-transformed fold change (log2FC), averaged expression across both sample groups, t values, p values, FDR-adjusted p values and differentially expressed classification of each gene in each sample pair.

1. **logfc**: Indicates the magnitude of change in gene expression. Negative values indicate down-regulation, while positive values indicate up-regulation.
2. **aveexpr**: The average expression level of the gene across all samples.
3. **t**: The t-statistic value used to assess the significance of the differential expression.
4. **p_value**: The raw p-value for the statistical test.
5. **adj_p_val**: The p-value adjusted for multiple comparisons.
6. **b**: The log-odds that the gene is differentially expressed.
7. **gene**: The Ensembl gene identifier.
8. **groupa**: The first group in the comparison.
  - Example: `D6` for the first comparison.
9. **groupb**: The second group in the comparison.
  - Example: `D5` for the second comparison.
10. **batch**: The batch identifier.
11. **compare**: The comparison being made (groupb/groupa).
  - Example: `D5/D6` indicates the comparison between D5 and D6.
12. **type**: Indicates whether the gene is up-regulated, down-regulated, or not differentially expressed (non-DEG).
  - Example: `down-regulate` indicates that the gene is down-regulated in the comparison.

#### the output of normality.R:
**normality_genelist_(date).csv**: The file contains p values ofv the Shapiro-Wilk test based on the ratio-based expression of each gene in the reference datasets.

1. **gene_compare**: Combination of gene identifier and comparison.
2. **log2fc_p**: P-value from the Shapiro-Wilk normality test for the log2FC values.
3. **gene**: Gene identifier.
4. **compare**: Comparison identifier (e.g., D5/D6).

#### the output of char_ratio.R:
**ref_expr_(date).csv**: The file contains the ratio-based expression of reference datasets, providing information on FC and median of limma-based p value for each gene.

1. **Var1**: Row identifier (index).
2. **Freq**: Frequency of detection across batches.
3. **gene**: Gene identifier.
4. **compare**: Comparison identifier (e.g., D5/D6).
5. **gene_compare**: Combined gene and comparison identifier.
6. **fc**: Fold change of the gene expression.
7. **medianp**: Median p-value for the gene comparison.

#### the output of char_refDEG.R:
**RefData_DEGs_(date).csv**: This file includes a list of genes classified as reference DEGs for each sample pair. Additionally, it lists the average FC and median p values across high-quality batches for each gene.

1. **var1**: Combined identifier of gene and comparison (e.g., "ENSG00000002919 D5/D6").
2. **freq**: Frequency of the gene appearing in the specified comparison.
3. **gene**: Gene identifier (e.g., "ENSG00000002919").
4. **compare**: Comparison group (e.g., "D5/D6").
5. **gene_compare**: Concatenation of gene identifier and comparison group for unique identification (e.g., "ENSG00000002919 D5/D6").
6. **fc**: Fold change value for the gene in the specified comparison.
7. **medianp**: Median p-value for the gene in the specified comparison, indicating the statistical significance of the differential expression.
8. **n_up**: Number of times the gene is classified as up-regulated in the comparison.
9. **n_non**: Number of times the gene is classified as non-DEG in the comparison.
10. **n_down**: Number of times the gene is classified as down-regulated in the comparison.
11. **final**: Final classification of the gene in the comparison. Possible values include:
    - **"non-DEG"**: Gene is not differentially expressed.
    - **"up-regulate"**: Gene is up-regulated.
    - **"down-regulate"**: Gene is down-regulated.
    - **"conflicting"**: Gene shows conflicting expression patterns (both up and down).

### Uncertainty
#### the output of uchar.R:
**Refdata_uchar_(date).csv**: This file contains the characterization uncertainty of each gene in the reference datasets.

1. **gene_compare**: Combined identifier of gene and comparison (e.g., "ENSG00000002919 D5/D6").
2. **gene**: Gene identifier (e.g., "ENSG00000002919").
3. **compare**: Comparison group (e.g., "D5/D6").
4. **uchar**: UCHAR value for the gene in the specified comparison. It represents the relative standard error of the fold change, multiplied by 100.

#### the output of ubb.R:
**Refdata_ubb_(date).csv**: This file contains the homogenity uncertainty of each gene in the reference datasets.

1. **gene**: The gene identifier (e.g., "ENSG00000001461").
2. **compare**: The comparison group, indicating the samples being compared (e.g., "F7/D6").
3. **gene_compare**: A concatenated identifier combining the gene and comparison group for unique identification (e.g., "ENSG00000001461 F7/D6").
4. **ubb**: The UBB value calculated for the given gene and comparison group, representing the unit-to-unit variability in the data.

#### the output of us.R:
**Refdata_us_(date).csv**: This file contains the instability uncertainty of each gene in the reference datasets.

1. **gene_compare**: A concatenated identifier combining the gene and comparison group for unique identification (e.g., "ENSG00000001461 F7/D6").
2. **gene**: The gene identifier (e.g., "ENSG00000001461").
3. **compare**: The comparison group, indicating the samples being compared (e.g., "F7/D6").
4. **us**: The US value calculated stability measure for the gene in the specified comparison.

#### the output of uc.R:
**RRefdata_uc_(date).csv**: This file contains the combined uncertainty of each gene in the reference datasets.

1. **Var1**: Combined identifier of gene and comparison (e.g., "ENSG00000002919 D5/D6").
2. **Freq**: Frequency of the gene appearing in the specified comparison.
3. **gene**: Gene identifier (e.g., "ENSG00000002919").
4. **compare**: Comparison group (e.g., "D5/D6").
5. **gene_compare**: Concatenation of gene identifier and comparison group for unique identification (e.g., "ENSG00000002919 D5/D6").
6. **fc**: Fold change value for the gene in the specified comparison.
7. **medianp**: Median p-value for the gene in the specified comparison, indicating the statistical significance of the differential expression.
8. **uchar**: Characterization uncertainty value.
9. **ubb**: Inhomogeneity uncertainty value.
10. **us**: Instability uncertainty value.
11. **uc**: Combined uncertainty value, calculated as the square root of the sum of the squares of uchar, ubb, and us.

#### the output of U.R:
**Refdata_U_(date).csv**: This file contains the extended uncertainty of each gene in the reference datasets.

1. **Var1**: Combined identifier of gene and comparison (e.g., "ENSG00000002919 D5/D6").
Freq: Frequency of the gene appearing in the specified comparison.
2. **gene**: Gene identifier (e.g., "ENSG00000002919").
3. **compare**: Comparison group (e.g., "D5/D6").
4. **gene_compare**: Concatenation of gene identifier and comparison group for unique identification (e.g., "ENSG00000002919 D5/D6").
5. **fc**: Fold change value for the gene in the specified comparison.
6. **medianp**: Median p-value for the gene in the specified comparison, indicating the statistical significance of the differential expression.
7. **uchar**: Characterization uncertainty value.
8. **ubb**: Inhomogeneity uncertainty value.
9. **us**: Instability uncertainty value.
10. **uc**: Combined uncertainty value, calculated as the square root of the sum of the squares of `uchar`, `ubb`, and `us`.
11. **U**: Extended uncertainty, calculated as `uc` multiplied by the coverage factor `k`.

### Performance metrics
#### the output of qc_rc.R:
**RC_each_batch_(date).csv**: This file contains the RC values between the test datasets and reference datasets for each batch.

1. **rc**: The correlation coefficient between the reference and test log fold changes for a specific batch. The values range from -1 to 1, where:

- 1 indicates a perfect positive correlation.
- 0 indicates no correlation.
- -1 indicates a perfect negative correlation.
2. **batch**: The identifier for the batch of samples used in the RNA-seq experiment. Each batch is labeled uniquely, indicating the platform and specific conditions used for sequencing (e.g., `P_BGI_L3_B1`).

#### the output of qc_mcc.R:
**MCC_each_batch_(date).csv**: This file contains the TP, TN, FP, FN, and MCC values between the test datasets and reference datasets for each batch.

1. **tp** (True Positives): The number of true positive predictions. These are cases where the model correctly predicts the positive class.
2. **tn** (True Negatives): The number of true negative predictions. These are cases where the model correctly predicts the negative class.
3. **fn** (False Negatives): The number of false negative predictions. These are cases where the model incorrectly predicts the negative class when the positive class is true.
4. **fp** (False Positives): The number of false positive predictions. These are cases where the model incorrectly predicts the positive class when the negative class is true.
5. **precision**: Precision is the ratio of true positive predictions to the total number of positive predictions (both true and false positives). 
6. **sensitivity** (Recall): Sensitivity, also known as recall, is the ratio of true positive predictions to the total number of actual positives (both true positives and false negatives). 
7. **specificity**: Specificity is the ratio of true negative predictions to the total number of actual negatives (both true negatives and false positives). 
8. **f1**: The F1 score is the harmonic mean of precision and recall. It provides a balance between precision and recall.
9. **mcc**: The Matthews correlation coefficient (MCC) is a measure of the quality of binary classifications. It takes into account true and false positives and negatives and is generally regarded as a balanced measure, even if the classes are of very different sizes.
10. **batch**: The identifier for the batch of samples used in the RNA-seq experiment. Each batch is labeled uniquely, indicating the platform and specific conditions used for sequencing (e.g., `P_BGI_L3_B1`).

#### the output of qc_rmse.R:
**RMSE_each_batch_(date).csv**: This file contains the RMSE value between test datasets and reference datasets for each batch.

1. **rmse**: The Root Mean Square Error (RMSE) value for the batch. This metric indicates the standard deviation of the residuals (prediction errors). A lower RMSE value indicates a better fit of the model to the reference data.

2. **batch**: The identifier for the batch of samples used in the RNA-seq experiment. Each batch is labeled uniquely, indicating the platform and specific conditions used for sequencing (e.g., P_BGI_L3_B1).


### Batch effects
#### the output of pvca.R:
(1) **PVCA_(date).csv**: This file contains the PVCA value of the expression profiles of multiple batches.

1. **name**:	The effect being analyzed (e.g., `batch`, `sample`)
2. **value**:	The weighted average proportion variance of the effect
3. **value2**:	Rounded value of the weighted average proportion variance
4. **type**:	Type of effect (`biological`, `technical`, `mixed`, `other`)

(2) **PVCA_plot_(date).pdf**: This file contains a bar plot of the PVCA results of the expression profiles of multiple batches.

#### the output of ratio.R:
**RNA_expression_ratio_(refsample)_(date).csv**: This file contains the ratio-based expression profiles of each gene.


## Reference

> Yu, Y., Hou, W., Liu, Y. et al. Quartet RNA reference materials improve the quality of transcriptomic data through ratio-based profiling. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01867-9
        
        
## Contact
If you have any questions, please feel free to contact us at quartet@fudan.edu.cn.

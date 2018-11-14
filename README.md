# rpd-rnaseq

### Introduction

**rpd-rnaseq** is a bioinformatics best-practice analysis pipeline for bulk RNAseq analyis and was
developed by the Research Pipeline Development Team at the [Genome Institute of
Singapore (GIS)](https://a-star.edu.sg/gis) based on a request by the Cancer Therapeutics & Stratified
Oncology 6.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across
multiple compute infrastructures in a very portable manner. It further comes with Conda, eocker and
Singularity support, making software installation easy and results reproducible. The implementation
tries to stay in line with [nf-core](https://nf-co.re/) practices and recommendations. The
pipeline itself was based on the [GIS STAR-RSEM pipeline
(2017-10)](https://github.com/gis-rpd/pipelines/blob/2017-10/rnaseq/star-rsem/README.md) and
[nf-core/rnaseq](https://github.com/nf-core/rnaseq) 



### Pipeline Steps

| Step                                                | Main program/s                      |
|-----------------------------------------------------|-------------------------------------|
| Trimmming, combining of read-pairs per sample and QC| Trim Galore, FastQC                 |
| Mapping                                             | STAR                                |
| Visualization files (bigWigs)                       | Samtools, deepTools (bamCoverage)   |
| RNASeq Metrics                                      | CollectRnaSeqMetrics (Picard)       |
| RNASeq QC                                           | RSeQC                               |
| Gene and isoform expression                         | RSEM                                |
| QC report                                           | MultiQC                             |


# Documentation


Documentation about the pipeline can be found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Usage](docs/usage.md)
3. [Output and results](docs/output.md)





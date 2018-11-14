# gis-rpd/rpd-rnaseq Output

## Overview

All results are written to the specified publish directory, which is `results` by default (which is
assumed below).

### MultiQC

Many results (FastQC, STAR alignment stats, RSeQC results) are summarised in MultiQC plots, which
can be found in `results/MultiQC/multiqc_report.html`

An [example MulltiQC reporte can be download here](multiqc_report.example.html.zip)

### Pipeline info

All pipeline internal info (timings etc.) are written to `results/pipeline_info/`. Most of this will not be
of interest to the average user. A good starting point is `results/pipeline_info/report.html`

### Result per sample

Results per sample can be found in the respectively named folders. If one of your samples is called
`ABC`, then you will find a folder called `results/ABC`, with all the following subfolders

#### bigwig

FIXME: describe

#### RnaSeqMetrics

See the [Picard RnaSeqMetric
Website](https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics) for
more information

FIXME: describe

#### RSEM

See the [RSEM website](https://github.com/deweylab/RSEM) for more information

FIXME: describe

#### RSeQC

RSeQC is a package of scripts designed to evaluate the quality of RNA seq data. You can find out more about the package at the [RSeQC website](http://rseqc.sourceforge.net/).

These are all quality metrics files and contains the raw data used for the plots in the MultiQC report. In general, the `.r` files are R scripts for generating the figures, the `.txt` are summary files, the `.xls` are data tables and the `.pdf` files are summary figures.

#### STAR

See [STAR website](https://github.com/alexdobin/STAR) for more information

* `Sample_Aligned.sortedByCoord.out.bam`
  * The aligned BAM file
* `Sample_Log.final.out`
  * The STAR alignment report, contains mapping results summary
* `Sample_Log.out` and `Sample_Log.progress.out`
  * STAR log files, containing a lot of detailed information about the run. Typically only useful for debugging purposes.
* `Sample_SJ.out.tab`
  * Filtered splice junctions detected in the mapping


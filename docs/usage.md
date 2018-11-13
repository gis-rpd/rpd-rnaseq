# rpd-rnaseq Usage

## General Nextflow info
Nextflow handles job submissions and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run gis-rpd/rpd-rnaseq -profile gis
```

This will launch the pipeline with the `gis` configuration profile.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

The work and results directory can be changed be changed on the commandline.

To get commandline usage informat, run
```bash

nextflow gis-rpd/rpd-rnaseq --help
```


### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull gis-rpd/rpd-rnaseq
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [gis-rpd/rpd-rnaseq releases page](https://github.com/gis-rpd/rpd-rnaseq/releases) and find the latest version number - numeric only (eg. `1.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main Arguments

To generate the online help, run

```bash
nextflow gis-rpd/rpd-rnaseq --help
```

### `-profile`

Use this parameter to choose a configuration profile. Each profile is designed for a different compute environment

* `gis`
    * Designed to be run on GIS' aquila
* `standard`
    * The default profile, used if `-profile` is not specified at all. Runs locally and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles.
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

### `-params-file`

This yaml files defines references and samples to be used.
See `params.test.yaml` for an example.

#### References section

This section lists all reference and indices being used by the pipeline.

```nextflow
references:
  staridx: /genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/
  star_gtf: /genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf
  refflat: /genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/refFlat.txt.gz
  gtf_bed: /genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.bed
  rsemidx: /genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/RSEMIndex/transcript
```

#### Samples section


This describes the FastQ input pair sample, broken up into so called read-units (think: read pairs). The second FastQ is optional.

The samples section is an exact copy of the yaml files which are create automatically for you in GIS when you are samples are demuxed. Note, currently two versions are circulating, one with readunits separated from samples (deprecated) and one with readunits listed under each samples. Use the latter.

The corresponding entries in params.yaml looks as follows:

```nextflow
samples:
  sample-name-1:
    readunits:
      readunit-1:
        fq1: path-to-R1.fastq.gz
        fq2: path-to-R2.fastq.gz
      ...
      readunit-n:
        fq1: path-to-R1.fastq.gz
        fq2: path-to-R2.fastq.gz
  ...
  sample-name-n:
    readunits:
      ..
}
```

So you can specify multiple samples and each samples can contain multiple fastq pairs (AKA readunits)

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:


### Library strandedness
Three command line flags / config parameters set the library strandedness for a run:

* `--forward_stranded`
* `--reverse_stranded`
* `--unstranded`

If not set, the pipeline will be run as unstranded. Specifying `--pico` makes the pipeline run in `forward_stranded` mode.

You can set a default in a cutom Nextflow configuration file such as one saved in `~/.nextflow/config` (see the [nextflow docs](https://www.nextflow.io/docs/latest/config.html) for more). For example:

```groovy
params {
    reverse_stranded = true
}
```

If you have a default strandedness set in your personal config file you can use `--unstranded` to overwrite it for a given run.

These flags affect the commands used in some stages of the pipeline

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.


### Adapter Trimming

To enable TrimGalore based trimming you have to enable it with `--trimming`
and then use any of the following command line parameters. These affect the command
used to launch TrimGalore

#### `--clip_r1 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).

#### `--clip_r2 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).

#### `--three_prime_clip_r1 [int]`
Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.

#### `--three_prime_clip_r2 [int]`
Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.


### Others

There are plenty more options that might be missing here. Let us know if anything is missing and also refer to the online usage:
```bash
nextflow gis-rpd/rpd-rnaseq --help
```


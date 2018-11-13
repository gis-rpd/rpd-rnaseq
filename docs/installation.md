# gis-rpd/rpd-rnaseq Installation

To start using the gis-rpd/rpd-rnaseq pipeline, there are three steps described below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)
3. Configure the pipeline

## 1) Install NextFlow

Nextflow is available on GIS' aquila.

Should you need to install it anyway, run the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if
`gis-rpd/rpd-rnaseq` is specified as the pipeline name.

## 3) Software dependencies

If you are within GIS and use `-profile GIS` then you don't have to do anything. For other cases,
you might have to prepare a Docker image or Singularity file or a conda environment (all
corresponding recipes are provided in this repo)

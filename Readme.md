# Snakemake workflow: Virome pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.14.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/rdenise/virome_pipeline/workflows/Tests/badge.svg?branch=main)](https://github.com/rdenise/virome_pipeline/actions?query=branch%3Amain+workflow%3ATests)
[![License (AGPL version 3)](https://img.shields.io/badge/license-GNU%20AGPL%20version%203-green.svg)](LICENSE)

## Aim

This software will produce ???

???

## Installation

### Step 1: install Snakemake and Snakedeploy

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

To install Mamba by conda, run

```shell
conda install mamba -n base -c conda-forge
```

Given that Mamba is installed, run 

```shell
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. 

#### Notes 

For all following commands (step 2 to 5) ensure that this environment is activated via 

```shell
conda activate snakemake
```

### Step 2: deploy workflow

 Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory on your system in the place of your choice as follow (note to change the path and file name to the one you want to create): : 

```shell
mkdir path/to/project-workdir
```

Then go to your project working directory as follow:

```shell
cd path/to/project-workdir
```

In all following steps, we will assume that you are inside of that directory.

Second, run 

```shell
snakedeploy deploy-workflow https://github.com/rdenise/virome_pipeline . --tag 1.0.0
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.

### Step 3: configure workflow

#### General settings

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.  

??? and ??? sheet
- ???

Missing values can be specified by empty columns or by writing `NA`.


### Step 4: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda, run Snakemake with 

```shell
snakemake --cores 1 --use-conda 
```

### Step 5: generate report

After finalizing your data analysis, you can automatically generate an interactive visual HTML report for inspection of results together with parameters and code inside of the browser via 

```shell
snakemake --report report.zip --report-stylesheet config/report.css
```
The resulting report.zip file can be passed on to collaborators, provided as a supplementary file in publications.

## Walk-Through and File Production

This pipeline consists of ??? steps called rules that take input files and create output files. Here is a description of the pipeline.

1. As snakemake is set up, there is a last rule, called `all`, that serves to call the last output files and make sure they were created.

2. A folder containing your work will be created:

```
   [project_name]/                           <- top-level project folder (your project_name)
   │
   │
   ├── report.html                           <- This report file      
   │
   ├── logs                                  <- Collection of log outputs, e.g. from cluster managers
   │
   ├── ????

```

3. Before starting the pipeline, your taxids will be checked and updated for the lowest descendant levels. This will create a file of checked taxids. It will skip this step if this file already exists.

4. When restarting the pipeline, the software will check if you made any changes in the seed file before running. If changes have been made, it will run what is necessary, else nothing will happen.

### Pipeline in image 

#### Normal behavior

<p align="center">
  <img src="doc/dummy_dag.png?raw=true" height="400">
</p>


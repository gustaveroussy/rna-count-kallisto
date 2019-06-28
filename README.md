Snakemake workflow: rna-count-kallisto

This is the [Snakemake](https://academic.oup.com/bioinformatics/article/28/19/2520/290322) workflow for RNA-Seq read count, powered by [Kallisto](https://pachterlab.github.io/kallisto/). Optional quality controls are performed by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregated by [MultiQC](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)

Each tool belong to their respective authors.

## 1. Install conda

In order to install conda, please follow the guide lines described on the official [Conda web-page](https://docs.continuum.io/anaconda/install/).

This will likely be:

```{sh}
# Conda pre-requisites
apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

# Download latest conda installer
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh

# Run the installer
bash ~/Downloads/Anaconda3-2019.03-Linux-x86_64.sh
```

## 2. Install Snakemake pre-requisites within Conda

We are going to use Python 3.7. No other options. Forget about python 2.7, it will be dropped out of support before 2020'. In [rna-count-kallisto](https://github.com/gustaveroussy/rna-count-kallisto), we use [formatted strings](https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting) available at python 3.6. It is not excluded that further development will include python 3.7 exclusive behaviours. Go for latest python !

However, old packages are still not available within [python package installer](https://pypi.org/). Especially, `datrie` fails to install within Python 3.7 pip module. So we need to install it first, then install Snakemake. That's also why there is not Conda environment for Snakemake with Python 3.7.

Here we go:

```{sh}
# Build virtual environment for Snakemake workflows
conda create -n SnakemakeWorkflows -c conda-forge datrie==0.7.1 python==3.7.3 pytest==4.5.0 git==2.21.0 drmaa==0.7.9

# Activate this environment
conda activate SnakemakeWorkflows

# Install Snakemake with the newly installed python
python3.7 -m pip install snakemake --user
```

Every single time you will use this pipeline, you have to activate that conda environment. The activation command will not be repeated everywhere in this wiki, so please, don't forget it!

## 3. Clone this workflow

The final step to use this workflow is to clone it to your computer. This can be easily done with git (which you have installed in your virtual environment few seconds ago!

Here we go:

```{sh}
git clone https://github.com/gustaveroussy/rna-count-kallisto.git
```

Everything is installed. Please, see the next page about [configurations](https://github.com/tdayris/rna-splicing-star-leafcutter/wiki/Configuration) before running. You can also see our [test] page if you want to test the [rna-splicing-star-leafcutter](https://github.com/tdayris/rna-splicing-star-leafcutter) pipeline.
## Local usage

Here is your command line:
```
snakemake
```

Yes. Nothing more as long as you hav [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/) and [Kallisto](https://pachterlab.github.io/kallisto/) in your `PATH`.

However, conda provides an efficient way to dynamically build virtual environments and to double check your tool versions, that's why we recommend:
```
snakemake --use-conda
```

And if you want to be sure that your OS is not interfering with your run, then use:
```
snakemake --use-conda --singularity
```

I won't detail more of the Snakemake possibilities, please see [their documentation](https://snakemake.readthedocs.io/en/stable/executable.html). However, many users like the following:
```
snakemake --use-conda --singularity --reason --printshellcmds
```

With these supplementary arguments, you will have the executed command lines printed to your STDOUT, and each rule will give you the reason why it is being executed. Nice, isn't it ?

## PBS Torque

In addition to the previous arguments, one might want to use this pipeline on a PBS Torque cluster.

```
snakemake --cluster " qsub -l walltime=00:{resources.time_min}:00,nodes=1:ppn={threads},mem={resources.mem_mb}mb -V " --jobname "{rulename}.snakejob.{jobid}" --use-conda --reason --printshellcmds --jobs 20  --restart-times 3
```

## Slurm

```
snakemake --use-conda --reason --printshellcmds --configfile config.yaml --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads} --time={resources.time_min}" --restart-times 3 --jobs 20
```
There are two configuration files: one for your analysis, one for describing your samples. Before each run, both of them are validated : their content is verified thanks to the [json-shcemas](https://json-schema.org/).

## 1. The general config file (yaml)

This is a [yaml](https://en.wikipedia.org/wiki/YAML) file. It could also be [json](https://en.wikipedia.org/wiki/JSON) formatted, however its name must not be changed. Yaml format has been preferred since
users requested a more human readable configuration file to write.

This configuration file contains the following sections (in any orders):

### ref

This is a very simple section: two values with evident significations. We need two paths, one to the fasta formatted genome sequence, and one other to the gtf formatted genome annotation. By GTF formatted, I mean GTF. No GFF, no GFF3, no GFF2, no GTF.gz, no GTF.bz2, nothing else than GTF. The same goes with the fasta file.

Example:
```
ref:
  fasta: path/to/fasta.fa
  gtf: /path/to/gtf.gtf
```

These paths can be either relative or absolute.

### params

This section contains additional command line arguments. One for all, **do not try to modify threading options** here. Do not try to modify  neither logging, nor temporary directories. These options are al
ready handled.

We have the following values:

copy_extra: Extra parameters for bash cp
kallisto_index_extra: Extra parameters for kallisto index rule
kallisto_quant_extra: Extra parameters for kallisto quant rule

Example:
```
params:
  copy_extra: "--parents --verbose"
  kallisto_index_extra: "--make-unique"
  kallisto_quant_extra: "--bias --bootstrap-samples=100"
```

### workflow

This section is used to activate or deactivate sections of the pipeline. In fact, if you just want to quantify and no quality control (you may have done them aside of the [rna-count-kallisto](https://github.com/gustaveroussy/rna-count-kallisto) pipeline), then turn these flag to `false` instead of `true`.

Warning, it is case sensitive!

Example:
```
workflow:
  fastqc: true
  multiqc: true
  aggregate: true
```

### General

The following parameters do not belong to any section, let them be at the top level of the yaml file:

design: The path to the design file
workdir: The path to the working directory
threads: The maximum number of threads used
singularity_docker_image: The image used within singularity
cold_storage: A list of cold storage mount points

Example:
```
design: design.tsv
workdir: .
threads: 1
singularity_docker_image: docker://continuumio/miniconda3:4.4.10
cold_storage:
  - /media
```

### Conclusion

A complete config.yaml file would look like this:

```
design: design.tsv
workdir: .
threads: 1
singularity_docker_image: docker://continuumio/miniconda3:4.4.10
cold_storage:
  - /media
ref:
  fasta: /path/to/genome/sequence.fa
  gtf: /path/to/genome/annotaton.gtf
workflow:
  fastqc: true
  multiqc: true
  aggregate: true
params:
  copy_extra: "--parents --verbose"
  kallisto_index_extra: "--make-unique"
  kallisto_quant_extra: "--bias --bootstrap-samples=100"
```

Remember, there is a python script here to build this file from command line!

## 2. The design file (tsv)

At this stage, we need a [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) file describing our analysis.

It must contain the following columns:

* Sample_id: the name of each samples
* Upstream_file: path to the upstream fastq file

The optional columns are:

* Downstream_file: path to downstream fastq files (usually your R2 in paired-end libraries)
* Any other information

Remember, there is a python script here to build this file from command line!
## Global workflow

This workflows takes fastq files, genome sequences and annotations as input, and returns abundance estimates along side with optional quality metrics.

## FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) stands for FastQ Quality Control. There is no publication associated with this tool, however, it remains an inescapable classic in bioinformatics.

This tool compiles a lot of informations about the raw reads obtained from sequencers. Actually, this tool does not perform any process **required** for the whole splicing analysis; however, it's a good practice.

Citation:
* Andrews, Simon. "FastQC: a quality control tool for high throughput sequence data." (2010).

## MultiQC

[MultiQC](https://multiqc.info/), just like FastQC, do not have any other purpose than quality metrics. It gathers all Flagstat and all FastQC individual metrics into one single report.

Citation:
* Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.

## Kallisto

[Kallisto](https://pachterlab.github.io/kallisto/) is a program for quantifying abundances of transcripts from RNA-Seq data. It is based on pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment. Pseudoalignment of reads preserves the key information needed for quantification, and kallisto is therefore not only fast, but also as accurate as existing quantification tools. In fact, because the pseudoalignment procedure is robust to errors in the reads, in many benchmarks kallisto significantly outperforms existing tools.

Citation:
* Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519

## Snakemake

[Snakemake](https://snakemake.readthedocs.io) is a pipeline/workflow manager written in python. It is used to handle the tools interaction, dependencies, command lines and cluster reservation. It is the skeleton of this pipeline. This pipeline is powered by the [Snakemake-Wrappers](https://snakemake-wrappers.readthedocs.io), the [Snakemake Workflows](https://github.com/snakemake-workflows), and the [conda](https://anaconda.org) project.

Citation:
* Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.

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

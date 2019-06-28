Snakemake workflow: rna-count-kallisto

This is the [Snakemake](https://academic.oup.com/bioinformatics/article/28/19/2520/290322) workflow for RNA-Seq read count, powered by [Kallisto](https://pachterlab.github.io/kallisto/). Optional quality controls are performed by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregated by [MultiQC](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)

Each tool belong to their respective authors.

Currently under active development, this pipeline is not intended to be used.

See the embedded [wiki](https://github.com/gustaveroussy/rna-count-kallisto/wiki) for more information.

For more detailed documentation, please see:

[Home](https://github.com/gustaveroussy/rna-count-kallisto/wiki)
* [Installation](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Installation)
  * [Conda](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Installation#1-install-conda)
  * [Snakemake](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Installation#2-install-snakemake-pre-requisites-within-conda)
  * [Workflow](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Installation#3-clone-this-workflow)
* [Configuration](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration)
  * [Config](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration#1-the-general-config-file-yaml)
    * [ref](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration#ref)
    * [params](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration#params)
    * [workflow](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration#workflow)
    * [general](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration#General)
  * [Desgin](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Configuration#2-the-design-file-tsv)
* [Usage](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Usage)
  * [Local](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Usage#local-usage)
  * [Torque](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Usage#pbs-torque)
  * [Slurm](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Usage#slurm)
* [Pipeline defails](https://github.com/gustaveroussy/rna-count-kallisto/wiki/Pipeline-details)

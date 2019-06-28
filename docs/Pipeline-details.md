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

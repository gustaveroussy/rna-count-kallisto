Material and Methods:
#####################

Quality control were made on raw `FastQ <https://en.wikipedia.org/wiki/FASTQ_format>`_ files with FastQC. Quality reports were gathered with MultiQC. Abundance estimation was performed with Kallisto. Optional parameters (if any) are listed below. Aggregation was performed with an in-house python script available in the `scripts` directory.

* Kallisto index optional arguments: `{{snakemake.config.params.kallisto_index_extra}}`
* Kallisto quantification optional arguments: `{{snakemake.config.params.kallisto_quant_extra}}`


Citations:
##########

MultiQC
  EWELS, Philip, MAGNUSSON, Måns, LUNDIN, Sverker, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016, vol. 32, no 19, p. 3047-3048.

  https://multiqc.info/

FastQC
  ANDREWS, Simon, et al. FastQC: a quality control tool for high throughput sequence data. 2010.

  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Kallisto
  Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519

  https://pachterlab.github.io/kallisto/

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

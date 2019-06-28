[rna-count-kallisto](https://github.com/gustaveroussy/rna-count-kallisto) is a [Snakemake Workflow](https://github.com/snakemake-workflows), it relies on the very same pre-requisites as [Snakemake](https://snakemake.readthedocs.io/) itself. Thus, it uses [pytest](https://docs.pytest.org/en/latest/) for its continuous testing. In addition, [Conda](https://anaconda.org/) is used to power all tools installations and virtual environment management.

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

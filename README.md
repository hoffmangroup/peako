# PeaKO

## What is peaKO?

PeaKO discovers motifs in ChIP-seq datasets with knockout controls. PeaKO takes in paired wild-type/knockout BAM files in addition to several reference files, as input. It returns a file of ranked motifs (see [our paper](https://doi.org/10.1093/bioadv/vbad031) for more details).


## Quick start

### Dependencies

1. Conda (Miniconda or Anaconda)
2. MEME Suite version 5.1.0 or MEME Suite version 4.12.0 with our CentriMo binary* (see below)

Please note that Conda can be installed locally without system administrator privileges. We suggest following [Miniconda's installation instructions](https://docs.conda.io/en/latest/miniconda.html).
PeaKO has only been tested on Linux systems.

### Installation

1. Download peaKO's [environment](https://peako.hoffmanlab.org/data/peako-env.yml) file.
2. Open a terminal and run `conda env create -f peako-env.yml` in your Downloads directory. This will create a Conda environment called "peako".
3. Run `conda activate peako` or `source activate peako` to activate this environment.
4. Install peaKO from PyPI by running `python3 -m pip install peako`.
5. You can test that this worked by running `peako --help`.

> **NOTES:** 
> If you run into errors indicating missing bs4 or pyYAML packages, try running `pip3 install beautifulsoup4==4.8.2 pyyaml`.
> If step 2 above fails to create a Conda environment, you may try:
>
>    conda create --name peako  
>    conda activate peako  # or source activate peako  
>    conda install python=3.7  
>    conda install -c anaconda beautifulsoup4=4.7 pandas  
>    conda install -c bioconda -c conda-forge -c anaconda snakemake-minimal flake8 pathlib2 ipython twine  
>    conda install -c bioconda pybedtools  
>
> Please note that we have only tested peaKO on Linux.

### Instructions for our modified CentriMo binary

If MEME Suite version 5.1.0 is installed on your system and accessible in your path, you do not need to install our CentriMo binary separately. 
Important: MEME Suite 5.1.0 must be installed from source and should not be installed from Conda at this time due to documented implementation issues.
Future versions of peaKO will support installation of MEME Suite through Conda once these issues are resolved.
If you are using an older version of MEME Suite, please follow the steps below to replace the CentriMo binary with our own to use peaKO.

1. Download MEME distribution 4.12.0 from the [MEME Suite Download page](http://meme-suite.org/doc/download.html).
2. Follow the "Quick Install" steps on the [MEME Suite Installation page](http://meme-suite.org/doc/install.html?man_type=web) up until `make install`.
4. After running `make install`, replace `$HOME/meme/bin/centrimo` with our [modified CentriMo binary](https://doi.org/10.5281/zenodo.3356995).
5. Make sure that `$HOME/meme/bin` is located on your `$PATH`. You should now be able to call `centrimo --help`.


## Usage

PeaKO uses Snakemake, which is a workflow management system.
You can run peaKO either locally or on a compute cluster using the Slurm job scheduling system.
To run on Slurm, you must create your own `cluster.config` file ([template](https://peako.hoffmanlab.org/data/cluster.json)) and provide it to peaKO via `--sm-cluster-config`.

Each step of the workflow either inherits from the main activated Conda environment ("peako") or uses its own separate environment.
If you are working on a compute cluster, run peaKO first with `--sm-build-envs` on a node with internet access to create these additional Conda environments.
Then, you can run it on the cluster without internet, providing a Slurm configuration file (see above).

After activating peaKO's Conda environment (`conda activate peako` or `source activate peako`), you can run peako as follows:

`peako <outdir> <wt-bam> <ko-bam> <organism> <chr-sizes> <trf-masked-genome> <motif-database> [options]`

There are 7 required arguments. Please provide absolute paths for files and directories.

- `outdir`: output directory (please make sure this already exists); all output directories and files will be created here
- `wt-bam`: wild-type sample BAM file
- `ko-bam`: knockout sample BAM file
- `organism`: name of organism (must be either `mouse` or `human`)
- `chr-sizes`: chromosome sizes file of reference genome (TXT)
- `trf-masked-genome`: TRF masked reference genome file (FASTA)
- `motif-database`: JASPAR motif database (MEME)

Here are the optional arguments:

General:

- `-h` or `--help`: access the help message and exit
- `-V` or `--version`: show the program's version and exit

PeaKO submodule:

- `-j <JASPAR_ID>`: transcription factor motif JASPAR identifier (e.g. MA0083.3)
- `-m <MOTIF>`: transcription factor motif common name (e.g. SRF)
- `--extra`: output all intermediate peaKO files for plotting
- `--pickle`: use pickled peaKO dictionaries from previous run 

Snakemake:

- `--sm-build-envs`: build conda environments for workflow and exit (requires internet connection)
- `--sm-cluster-config`: snakemake cluster configuration file (JSON)


## Output

Currently, peaKO generates output directories and files for each step.
These can all be found under your provided `outdir` directory.
PeaKO's main output file is `<outdir>/peako_out/peaKO-rankings.txt`, which contains a ranked list of motifs.


## Additional resources

Source code is available at: https://github.com/hoffmangroup/peako.

We have deposited the [current version of the code](https://doi.org/10.5281/zenodo.3338330), [example HTML and TXT CentriMo outputs](https://doi.org/10.5281/zenodo.3338324), and a [modified CentriMo binary](https://doi.org/10.5281/zenodo.3356995) on Zenodo.


## Citation

If you found peaKO useful, please cite:



[Danielle Denisko *, Coby Viner *, Michael M Hoffman, Motif elucidation in ChIP-seq datasets with a knockout control, *Bioinformatics Advances* 2023, vbad031](https://doi.org/10.1093/bioadv/vbad031).

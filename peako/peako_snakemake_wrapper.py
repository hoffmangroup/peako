#!/usr/bin/env python

import argparse
import os
import logging
import sys
import subprocess
import json

from pkg_resources import resource_filename
#from .version import __version__

snakefile = resource_filename(__name__, "data/snakefile")


def get_basename_without_exts(file_path):
    # adapted from: https://stackoverflow.com/questions/678236/
    file_basename = os.path.basename(file_path)
    filename_without_extension = file_basename.split('.')[0]
    return filename_without_extension


def get_sample_name(sample):
    if sample.endswith(('.bam.gz', '.bam')):
        if sample.endswith('.bam'):
            sample_name = get_basename_without_exts(sample)
    else:
        sys.exit("{} must end in .bam or .bam.gz".format(sample))
    return sample_name


def check_and_create_symlink(directory, original_filepath, symlink_filepath):
    if not os.path.isdir(directory):
        subprocess.run("mkdir {}".format(directory), shell=True)
    if not os.path.isfile(symlink_filepath):
        subprocess.run("ln -s {} {}".format(original_filepath,
                                            symlink_filepath),
                       shell=True)
    else:
        logging.info("{} already exists".format(symlink_filepath))


def main():
    logging.basicConfig(format='%(asctime)s >peaKO-workflow: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog='peako',
        description="""
PeaKO ranks motifs from ChIP-seq experiments performed with knockout
controls. This peaKO workflow processes a pair of wild-type/knockout
BAM files and returns motif rankings. PeaKO computes its ranking
metric by comparing CentriMo results across two different ways of
differentially accounting for the knockout dataset (Pipeline A vs.
Pipeline B; see manuscript for details).
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

Example CentriMo HTML and TXT files are available at: \
https://doi.org/10.5281/zenodo.3338330.

Usage example:
$ peako ~/workdir/peako-workflow/ \

        ~/analysis/tead4-wt.bam   \

        ~/analysis/tead4-ko.bam   \

        mouse                     \

        ~/mm10/chrsize.txt         \

        ~/trf/mm10-masked.fa       \

        ~/motif_database/JASPAR_CORE_2016_vertebrates.meme

If you found peaKO useful, please cite:

Denisko D, Viner C, Hoffman MM. Motif elucidation in ChIP-seq datasets with a
knockout control. BioRxiv <ID> [Preprint]. 2019. \
Available from: https://doi.org/<ID>
        """)

    parser.add_argument('workdir', metavar="outdir",
                        help='output directory')
    parser.add_argument('wt_bam', metavar="wt-bam",
                        help='wild-type sample BAM file')
    parser.add_argument('ko_bam', metavar="ko-bam",
                        help='knockout sample BAM file')

    parser.add_argument('organism', metavar="organism",
                        help='name of organism (mouse or human)')
    parser.add_argument('chr_sizes', metavar="chr-sizes",
                        help='chromosome sizes file of reference genome (TXT)')
    parser.add_argument('trf_genome', metavar="trf-masked-genome",
                        help='TRF masked reference genome file (FASTA)')
    parser.add_argument('motif_database', metavar="motif-database",
                        help='motif database (MEME)')

    # general
#    parser.add_argument('-V', '--version', action='version',
#                        version=__version__)

    # peaKO submodule
    parser.add_argument('-j', dest='jaspar_id',
                        help='transcription factor motif JASPAR ' +
                        'identifier (e.g. MA0083.3)')
    parser.add_argument('-m', dest='motif',
                        help='transcription factor motif common name ' +
                        '(e.g. SRF)')
    parser.add_argument('--extra', action="store_true",
                        help='output all intermediate peaKO files for ' +
                        'plotting')
    parser.add_argument('--pickle', action="store_true",
                        help='use pickled peaKO dictionaries from previous ' +
                        'run')

    # snakemake
    parser.add_argument('--sm-build-envs', action="store_true",
                        help='build conda environments for workflow and ' +
                        'exit (requires internet connection)')
    parser.add_argument('--sm-cluster-config', dest='sm_cluster_config',
                        help='snakemake cluster configuration file (JSON)')
    # parser.add_argument('--sm-config', dest='sm_config',
    #                     help='snakemake configuration file (JSON)')

    args = parser.parse_args()
    # logging.info('Command line arguments: {}'.format(args))

    workdir = args.workdir
    wt_bam = args.wt_bam
    ko_bam = args.ko_bam
    organism = args.organism
    motif_database = args.motif_database
    chr_sizes = args.chr_sizes
    trf_genome = args.trf_genome
    peako_jaspar_id = args.jaspar_id
    peako_motif = args.motif
    peako_extra = args.extra
    peako_pickle = args.pickle

    script_path = os.path.dirname(os.path.realpath(__file__))

    # check that directory and all files exist
    if not os.path.isdir(workdir):
        sys.exit("{} is not a directory... exiting!".format(workdir))

    file_list = [wt_bam, ko_bam, chr_sizes, trf_genome]
    for f in file_list:
        if not os.path.isfile(f):
            sys.exit("{} is not a file... exiting!".format(f))

    # derive sample names
    wt_sample = get_sample_name(wt_bam)
    ko_sample = get_sample_name(ko_bam)
    pipb_sample = "-".join([wt_sample, ko_sample])

    # create symlinks in workdir/data and workdir/reference_files
    # (if not already present)
    data_dir = os.path.join(workdir, "data")
    ref_dir = os.path.join(workdir, "reference_files")

    wt_symlink = os.path.join(data_dir, "{}-wt.bam".format(wt_sample))
    ko_symlink = os.path.join(data_dir, "{}-ko.bam".format(ko_sample))
    motif_db_symlink = os.path.join(ref_dir, os.path.basename(motif_database))
    chr_sizes_symlink = os.path.join(ref_dir, os.path.basename(chr_sizes))
    trf_genome_symlink = os.path.join(ref_dir, os.path.basename(trf_genome))

    check_and_create_symlink(data_dir, wt_bam, wt_symlink)
    check_and_create_symlink(data_dir, ko_bam, ko_symlink)
    check_and_create_symlink(ref_dir, motif_database, motif_db_symlink)
    check_and_create_symlink(ref_dir, chr_sizes, chr_sizes_symlink)
    check_and_create_symlink(ref_dir, trf_genome, trf_genome_symlink)

    # associate organism with mm or hs for MACS2
    if organism == "mouse":
        macs2_genome_size = "mm"
    elif organism == "human":
        macs2_genome_size = "hs"
    else:
        sys.exit("{} is not a valid organism for ".format(organism) +
                 "this workflow (must choose mouse or human)... exiting!")

    # convert peako args to a string
    peako_params = []
    if peako_extra:
        peako_params.append("--extra")
    if peako_pickle:
        peako_params.append("--pickle")
    if peako_jaspar_id:
        peako_params.append("-j {}".format(peako_jaspar_id))
    if peako_motif:
        peako_params.append("-m {}".format(peako_motif))

    if peako_params:
        peako_params_command = " ".join(peako_params)
    else:
        peako_params_command = ""
    
    # convert snakemake args to string
    # sm_additional_params = []
    # if sm_unlock:
    #     sm_additional_params.append("--unlock")

    # write everything to JSON file
    config_file = os.path.join(workdir, "config.json")
    with open(config_file, "w") as f:
        json.dump({'sample_wt_prefix': wt_sample,
                   'sample_ko_prefix': ko_sample,
                   'sample_pipb_prefix': pipb_sample,
                   'macs2_genome_size': macs2_genome_size,
                   'macs2_q_value': 0.05,
                   'chr_sizes': ("reference_files/{}"
                                 .format(os.path.basename(chr_sizes))),
                   'trf_masked_genome': ("reference_files/{}"
                                         .format(os.path
                                                 .basename(trf_genome))),
                   'motif_database': ("reference_files/{}"
                                      .format(os.path
                                              .basename(motif_database))),
                   'min_meme_motif_width': 7,
                   'max_meme_motif_width': 12,
                   'peako_params': peako_params_command,
                   'script_path': script_path},
                  f, indent=4)

    # build conda environments and exit
    if args.sm_build_envs:
        subprocess.run('snakemake -s {} -d {} '.format(snakefile, workdir) +
                       '--use-conda --create-envs-only',
                       shell=True)
        sys.exit("Built Conda environments for peaKO. You can now run " +
                 "it on a cluster without internet if you wish.")

    # source activate conda enrironment (has to already be built)
    if args.sm_cluster_config:
        cluster_command = ('"sbatch -p {cluster.partition} ' +
                           '-t {cluster.time} --mem {cluster.mem} ' +
                           '-J {cluster.name} -o {cluster.log}"')
        subprocess.run('snakemake -s {} -d {} '.format(snakefile, workdir) +
                       '--configfile {} '.format(config_file) +
                       '-j 4 --latency-wait 60 --use-conda ' +
                       '--cluster-config {} '.format(args.sm_cluster_config) +
                       '--cluster {}'.format(cluster_command),
                       shell=True)
    else:
        subprocess.run('snakemake -s {} -d {} '.format(snakefile, workdir) +
                       '--configfile {} --use-conda '.format(config_file) +
                       '--latency-wait 60',
                       shell=True)


if __name__ == '__main__':
    main()

# TO DO:
# change verbosity of all programs at once
# maybe give option to provide custom config and skip over all of these steps

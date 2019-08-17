#!/usr/bin/env python

"""
PeaKO ranks motifs based on peak overlap between differential
MEME-ChIP (DMC) and knockout implemented normalization (KOIN)
pipeline peak sets.

Note: In our manuscript, we apply the following mappings:
     DMC pipeline --> Pipeline A
     KOIN pipeline --> Pipeline B
"""

from __future__ import absolute_import, division, print_function
import argparse
import linecache
from math import log
from collections import OrderedDict
import csv
import os
from itertools import chain
import logging

import yaml
from bs4 import BeautifulSoup
import pandas as pd
import pybedtools
import pickle
from six.moves import range


class CentriMoOutput:

    def __init__(self, filename):
        self.filename = filename

    def data(self):
        return store_data_from_html(self.filename)

    def num_db_motifs(self):
        """
        Number of motifs in provided databases
        (usually JASPAR, MEME, and/or DREME).
        """
        return sum([db['count'] for db in self.data()['motif_dbs']])

    def num_motifs(self):
        """
        Number of motifs that CentriMo found to be sufficiently
        enriched in peak sets.
        """
        return len(self.data()['motifs'])


def store_data_from_html(filename_html):
    """
    Return dictionary of HTML file between script tags.

    Use Beautiful Soup module from Python to extract text between script tags
    in HTML file. Then, split on '=' and remove semicolons.

    We are extracting the set of all motifs reported by CentriMo.


    Keyword arguments:
    filename_html -- full path of HTML file
    """
    with open(filename_html, 'r') as f:
        html = f.read()
        soup = BeautifulSoup(html, 'html.parser')
        script = soup.find('script')
        # adapted from: http://stackoverflow.com/questions/13323976/
        json_text = script.text.split('=', 1)[1].strip().rstrip(';')
        data = yaml.safe_load(json_text)  # loads into string objects
        return data


def read_centrimo_df(filename_txt):
    """
    Read CentriMo TXT file and return dataframe.
    """
    line = linecache.getline(filename_txt, 2)  # line 1: warning
    header = line.split()[1:]  # split on whitespace and ignore leading '#'
    df = pd.read_csv(filename_txt, header=None, names=header,
                     skip_blank_lines=True, comment='#', delim_whitespace=True)
    return df


def read_and_sort_centrimo_df(filename_txt):
    """
    Read CentriMo TXT file and return dataframe sorted by (log-adjusted)
    binomial p-values (most significant motifs at top).


    Keyword arguments:
    filename_txt -- full path of CentriMo TXT file
    """
    df = read_centrimo_df(filename_txt)
    sorted_df = (df
                 .sort_values(by='log_adj_p-value', ascending=True)
                 .reset_index(drop=True))
    sorted_df.index += 1
    return sorted_df


def read_and_sort_fisher_exact_centrimo_df(filename_txt):  # used in counts
    """
    Read differential CentriMo TXT file and return dataframe sorted by
    (log-adjusted) Fisher's Exact test p-values (most significant motifs
    at top).

    Key word arguments:
    filename_txt -- full path of differential CentriMo TXT file
    """
    df = read_centrimo_df(filename_txt)
    sorted_df = (df
                 .sort_values(by='log_fisher_adj_pvalue', ascending=True)
                 .reset_index(drop=True))
    sorted_df.index += 1
    return sorted_df


def filter_centrimo_df(df, CentriMoOutput_):
    """
    Filter out motifs with E-values greater than 0.1.


    Keyword arguments:
    df -- sorted dataframe of CentriMo TXT values
    CentriMoOutput_ -- CentriMoOutput class object
    """
    log_adj_p_val_filter = log(0.1/CentriMoOutput_.num_db_motifs())
    filtered_motif_list = (df
                           .loc[df['log_adj_p-value'] <=
                                log_adj_p_val_filter]['id'])
    return filtered_motif_list


def write_dict_to_bed_file(dictionary, outdir, file_list, filetype):
    """
    Write contents of peak set dictionary to BED file.


    Keyword arguments:
    dictionary -- keys are motifs, values are genomic region intervals
    file_list -- list to append newly created BED files
    filetype -- one of 'intersection', 'koin', 'neg', 'pos'
    """
    for motif, regions in dictionary.items():

        if filetype == 'intersection':
            bed_file = os.path.join(outdir, 'intcount-{}.bed'.format(motif))
            file_list.append(bed_file)
            with open(bed_file, 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter='\t',
                                    lineterminator='\n')
                for reg_index in range(0, len(regions)):
                    region = regions[reg_index]
                    writer.writerow([region[0], region[1], region[2]])
        else:
            if filetype == 'koin':
                bed_file = os.path.join(outdir, 'koin-{}.bed'.format(motif))
            elif filetype == 'neg':
                bed_file = os.path.join(outdir, '{}-neg.bed'.format(motif))
            elif filetype == 'pos':
                bed_file = os.path.join(outdir, '{}.bed'.format(motif))

            file_list.append(bed_file)
            with open(bed_file, 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter='\t',
                                    lineterminator='\n')
                for reg_index in range(0, len(regions)):
                    chrm, region = regions[reg_index].split(':')
                    start, end = region.split('-')
                    writer.writerow([chrm, int(start), int(end)])
        logging.debug('Wrote {} file.'
                      .format(os.path.basename(bed_file)))


def pickle_dictionary(dictionary, pickle_file):
    with open(pickle_file, 'wb') as handle:
        pickle.dump(dictionary, handle)


def load_pickled_dict(pickle_file):
    with open(pickle_file, 'rb') as handle:
        return pickle.load(handle)


def remove_prefix(text, prefix):
    # adapted from https://stackoverflow.com/questions/16891340
    if text.startswith(prefix):
        return text[len(prefix):]


def main():
    logging.basicConfig(format='%(asctime)s >peaKO: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog='peako',
        description="""
        PeaKO ranks motifs from ChIP-seq experiments performed with knockout
        controls. PeaKO computes its ranking metric by comparing CentriMo
        results across two different ways of differentially accounting for
        the knockout dataset (Pipeline A vs. Pipeline B; see manuscript
        for details).

        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

Example CentriMo HTML and TXT files are available at: \
https://doi.org/10.5281/zenodo.3338330.

Usage example:
$ ./peako.py ~/workdir/centrimo/centrimo-a/ \

             centrimo.html                  \

             centrimo.txt                   \

             ~/workdir/centrimo/centrimo-b/ \

             centrimo.html                  \

             centrimo.txt                   \

             ~/workdir/peako/

If you found peaKO useful, please cite:

Denisko D, Viner C, Hoffman MM. Motif elucidation in ChIP-seq datasets with a
knockout control. BioRxiv <ID> [Preprint]. 2019. \
Available from: https://doi.org/<ID>
        """)
    # req_grp = parser.add_argument_group(title='required arguments')
    parser.add_argument('indir_a', metavar='pipeline-a-dir',
                        help='differential CentriMo directory ' +
                        '(Pipeline A)')
    parser.add_argument('html_a', metavar='pipeline-a-html',
                        help='differential CentriMo HTML filename ' +
                        '(Pipeline A)')
    parser.add_argument('txt_a', metavar='pipeline-a-txt',
                        help='differential CentriMo TXT filename ' +
                        '(Pipeline A)')
    parser.add_argument('indir_b', metavar='pipeline-b-dir',
                        help='non-differential CentriMo directory ' +
                        '(Pipeline B)')  # KOIN
    parser.add_argument('html_b', metavar='pipeline-b-html',
                        help='non-differential CentriMo HTML filename ' +
                        '(Pipeline B)')  # KOIN
    parser.add_argument('txt_b', metavar='pipeline-b-txt',
                        help='non-differential CentriMo TXT filename ' +
                        '(Pipeline B)')  # KOIN
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0')
    parser.add_argument('-v', '--verbose', action="store_true",
                        help='print detailed progress')
    parser.add_argument('-m', dest='motif',
                        help='transcription factor motif common name ' +
                        '(e.g. SRF)')
    parser.add_argument('-j', dest='jaspar_id',
                        help='transcription factor motif JASPAR ' +
                        'identifier (e.g. MA0083.3)')
    parser.add_argument('--pickled', action="store_true",
                        help='use pickled dictionaries from previous run(s) ' +
                        'located in the output directory (outdir)')
    parser.add_argument('--extra', action="store_true",
                        help='output all intermediate files for plotting')
    parser.add_argument('--filter', action="store_true",
                        help='filter out half of peaKO motifs, according to ' +
                        'the number of intersecting regions')

    args = parser.parse_args()
    logging.info('Command line arguments: {}'.format(args))

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # dmc_html = os.path.join(args.d, args.d1)
    dmc_html = os.path.join(args.indir_a, args.html_a)
    # dmc_txt = os.path.join(args.d, args.d2)
    dmc_txt = os.path.join(args.indir_a, args.txt_a)
    # koin_html = os.path.join(args.k, args.k1)
    koin_html = os.path.join(args.indir_b, args.html_b)
    # koin_txt = os.path.join(args.k, args.k2)
    koin_txt = os.path.join(args.indir_b, args.txt_b)
    outdir = args.outdir
    transcription_factor = args.motif
    jaspar_id = args.jaspar_id
    pickled = args.pickled
    extra = args.extra

    # Create CentriMoOutput type objects to store HTML.
    dmc_html_obj = CentriMoOutput(dmc_html)
    koin_html_obj = CentriMoOutput(koin_html)
    logging.info('Read HTML files.')

    # Create dataframes from TXT files and sort according to log-adjusted
    # (binomial) p-values.
    sorted_dmc_df = read_and_sort_centrimo_df(dmc_txt)
    sorted_koin_df = read_and_sort_centrimo_df(koin_txt)

    # Filter out motifs with E value > 0.1 and
    # append those remaining to a list.
    filtered_dmc_motif_list = filter_centrimo_df(sorted_dmc_df,
                                                 dmc_html_obj).tolist()
    filtered_koin_motif_list = filter_centrimo_df(sorted_koin_df,
                                                  koin_html_obj).tolist()
    logging.info('Created and sorted dataframes. Filtered out motifs'
                 ' with E > 0.1. Pipeline A list contains {} motifs.'
                 ' Pipeline B list contains {} motifs.'.format(
                     len(filtered_dmc_motif_list),
                     len(filtered_koin_motif_list)))

    # Create dictionaries:
    # keys: motifs passing E-value filtering
    # values: chromosome regions where these motifs bind
    #         (as determined by CentriMo)
    #         (format: lists containing 'chri:xxxxxxx-xxxxx5xx' strings)

    dct = OrderedDict()
    neg_dct = OrderedDict()
    koin_dct = OrderedDict()

    # Calculate total number of motifs and make lists of motif IDs from
    # original HTML files.
    number_of_motifs = dmc_html_obj.num_motifs()
    dmc_pos_html_motifs = [d['id'] for d in dmc_html_obj.data()['motifs']]
    number_of_koin_motifs = koin_html_obj.num_motifs()
    koin_html_motifs = [d['id'] for d in koin_html_obj.data()['motifs']]
    logging.info('Calculated total number of motifs and created ID lists.')

    # Make dictionaries to store motifs or load in pickled dicts
    dct_pickle_file = os.path.join(outdir, 'dct.pickle')
    neg_dct_pickle_file = os.path.join(outdir, 'neg_dct.pickle')
    koin_dct_pickle_file = os.path.join(outdir, 'koin_dct.pickle')

    if not pickled:
        # Create dictionaries.
        logging.warning('No pickled dictionaries, so creating some.')

        # Make lists of all sequences from HTMLs.
        pos_sequence_list = dmc_html_obj.data()['sequences']
        neg_sequence_list = dmc_html_obj.data()['neg_sequences']
        koin_sequence_list = koin_html_obj.data()['sequences']
        logging.info('Created sequence lists.')

        for filtered_dmc_motif in filtered_dmc_motif_list:

            for index, motif in enumerate(dmc_pos_html_motifs):

                if motif == filtered_dmc_motif:
                    logging.debug('Creating dict entry for {}...'.format(
                        motif))
                    # array of integers (indices) correspond
                    # to regions in dmc_html_obj.data()['sequences']
                    motif_dct = dmc_html_obj.data()['motifs'][index]

                    for seqnum in motif_dct['seqs']:
                        if "id{}".format(filtered_dmc_motif) in dct:
                            # adapted from http://stackoverflow.com/
                            # questions/3199171
                            dct["id{}"
                                .format(filtered_dmc_motif)].append(
                                    pos_sequence_list[seqnum])
                        else:
                            dct["id{}"
                                .format(filtered_dmc_motif)] = [
                                    pos_sequence_list[seqnum]]

                    for neg_seqnum in motif_dct['neg_seqs']:
                        if "id{}".format(filtered_dmc_motif) in neg_dct:
                            neg_dct["id{}".format(filtered_dmc_motif)].append(
                                neg_sequence_list[neg_seqnum])
                        else:
                            neg_dct["id{}".format(filtered_dmc_motif)] = [
                                neg_sequence_list[neg_seqnum]]
                        break
        logging.info('Created Pipeline A positive and negative sequence ' +
                     'dictionaries...')

        for filtered_koin_motif in filtered_koin_motif_list:

            for index, motif in enumerate(koin_html_motifs):

                if motif == filtered_koin_motif:
                    logging.debug('Creating Pipeline B dict entry ' +
                                  'for {}...'.format(motif))
                    motif_dct = koin_html_obj.data()['motifs'][index]

                    for seqnum in motif_dct['seqs']:
                        if "id{}".format(filtered_koin_motif) in koin_dct:
                            koin_dct["id{}".format(
                                filtered_koin_motif)].append(
                                    koin_sequence_list[seqnum])
                        else:
                            koin_dct["id{}".format(filtered_koin_motif)] = [
                                koin_sequence_list[seqnum]]
                    break
        logging.info('... and Pipeline B sequence dictionaries...')
        logging.info('Created dictionaries to store motifs and '
                     'their associated chromosome regions.')

        # Create 3-column BED files for motifs in each dictionary.
        pickle_dictionary(dct, dct_pickle_file)
        pickle_dictionary(neg_dct, neg_dct_pickle_file)
        pickle_dictionary(koin_dct, koin_dct_pickle_file)

        logging.info('Pickled dictionaries. ')

    elif pickled:
        logging.info('Loading pickled dictionaries...')
        # to reload pickled objects
        dct = load_pickled_dict(dct_pickle_file)
        neg_dct = load_pickled_dict(neg_dct_pickle_file)
        koin_dct = load_pickled_dict(koin_dct_pickle_file)
        logging.info('Loaded!')
    # else:
        # logging.warning('Neither pickled nor not pickled...')

    dmc_file_list = list()
    negMCfile_list = list()
    koinfile_list = list()

    write_dict_to_bed_file(dct, outdir, dmc_file_list, 'pos')
    write_dict_to_bed_file(neg_dct, outdir, negMCfile_list, 'neg')
    write_dict_to_bed_file(koin_dct, outdir, koinfile_list, 'koin')

    logging.info('Wrote dictionaries to BED files. ')

    # Calculate ratio score and print results to file.

    intersection_regs_dct = OrderedDict()
    peako_file_unsorted = os.path.join(outdir, 'peaKO-rankings-unsorted.txt')

    with open(peako_file_unsorted, 'w') as csv_file:
        writer = csv.writer(csv_file, delimiter='\t', lineterminator='\n')
        # number_of_intersection_regs, number_of_DMC_positive_regs,
        # number_of_K_regs
        writer.writerow(["motif_ID", "intersection_regions",
                         "pipeline_a_wt_regions", "pipeline_b_regions",
                         "ratio", "alt_ID"])

        for dmc_file in dmc_file_list:
            dmc_pos_bedtool = pybedtools.BedTool(dmc_file)
#            motif_id = os.path.basename(dmc_file)[:-4]  # assume *.bed # XXX
            dmc_base, dmc_ext = os.path.splitext(dmc_file)
            assert dmc_ext == ".bed"
            motif_id = os.path.basename(dmc_base)
            dmc_neg_file = os.path.join(outdir, '{}-neg.bed'.format(motif_id))

            for koin_file in koinfile_list:
                koin_base, koin_ext = os.path.splitext(koin_file)
                assert koin_ext == ".bed"
                if remove_prefix(os.path.basename(koin_base),
                                 'koin-') == motif_id:
                    # if os.path.basename(koin_file)[5:-4] == motif_id:
                    koin_bedtool = pybedtools.BedTool(koin_file)

                    try:
                        dmc_neg_bedtool = pybedtools.BedTool(dmc_neg_file)

                        # 'A=True' removes full feature in dmc_pos_bedtool
                        # if overlapping
                        dmc_difference = dmc_pos_bedtool.subtract(
                            dmc_neg_bedtool, A=True)

                        # 'wa=True' reports full feature of dmc_difference
                        # if overlapping
                        koin_intersect_dmc = (dmc_difference
                                              .intersect(koin_bedtool,
                                                         wa=True))
                    # workaround for no negative file
                    except (FileNotFoundError, ValueError):
                        logging.warning('The negative sequences bed file,'
                                        ' {}, does not exist.'
                                        .format(os.path.basename(
                                            dmc_neg_file)))
                        koin_intersect_dmc = (dmc_pos_bedtool
                                              .intersect(koin_bedtool,
                                                         wa=True))

                    koin_intersect_dmc_merged = (koin_intersect_dmc
                                                 .sort()
                                                 .merge())
                    number_of_intersection_regs = (koin_intersect_dmc_merged
                                                   .count())
                    logging.debug('Motif {} has {} intersection ' +
                                  'regions'.format(
                                      motif_id, number_of_intersection_regs))

                    # Create dictionary entry to print intersecting regions
                    # for each motif.
                    for interval in koin_intersect_dmc_merged:
                        if "{}".format(motif_id) in intersection_regs_dct:
                            (intersection_regs_dct["{}".format(motif_id)]
                             .append(interval))
                        else:
                            intersection_regs_dct["{}"
                                                  .format(
                                                      motif_id)] = [interval]

                    # Find alternate ID using original HTML file.
                    for motif in dmc_html_obj.data()['motifs']:
                        if motif_id == ('id' + motif['id']):
                            try:
                                motif_id_alt = motif['alt']
                            except (FileNotFoundError, ValueError):
                                motif_id_alt = motif['id']

                            ratio = float(
                                number_of_intersection_regs)/len(dct[motif_id])
                            writer.writerow([motif_id,
                                             number_of_intersection_regs,
                                             len(dct[motif_id]),
                                             len(koin_dct[motif_id]),
                                             ratio, motif_id_alt])

    if not extra:
        for dmc_file in dmc_file_list:
            os.remove(dmc_file)
            dmc_basename, dmc_ext = os.path.splitext(
                os.path.basename(dmc_file))
            dmc_neg_file = os.path.join(os.path.dirname(dmc_file),
                                        dmc_basename + '-neg' + dmc_ext)
            if os.path.isfile(dmc_neg_file):
                os.remove(dmc_neg_file)
        for koin_file in koinfile_list:
            os.remove(koin_file)
    logging.info('Wrote first (unsorted and unfiltered) ratio file.')

    # Print out intersection regions.
    if extra:
        intersection_regs_list = list()
        write_dict_to_bed_file(intersection_regs_dct, outdir,
                               intersection_regs_list, 'intersection')
        logging.info('Wrote BED files of intersection regions.')

        intersection_counts_list = (set(list(chain
                                             .from_iterable(
                                                 list(intersection_regs_dct.values())))))
        intersection_counts_file = os.path.join(outdir,
                                                'all-intcount-regions.bed')

        with open(intersection_counts_file, 'w') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t', lineterminator='\n')
            for region in intersection_counts_list:
                writer.writerow([region[0], region[1], region[2]])
        logging.info('Wrote all intersection regions to file.')

    # Filter out the bottom 50% of motifs when ranked by the number
    # of intersecting regions.
    ratio_filepath = peako_file_unsorted
    rankings_df = pd.read_csv(ratio_filepath, delim_whitespace=True)
    os.remove(peako_file_unsorted)
    rankings_ratio_sorted_df = (rankings_df
                                .sort_values(by=['ratio',
                                                 'intersection_regions',
                                                 'pipeline_a_wt_regions',
                                                 'pipeline_b_regions'],
                                             ascending=False)
                                .reset_index(drop=True))
    rankings_ratio_sorted_df.index += 1

    peako_file = os.path.join(outdir, 'peaKO-rankings.txt')
    rankings_ratio_sorted_df.to_csv(peako_file,
                                    index=True, sep='\t')
    logging.info('Wrote peaKO file (sorted).')

    filtering_threshold = int(len(rankings_df)/2)  # 50% of motifs
    rankings_intcount_filtered_df = (rankings_df
                                     .sort_values(by=['intersection_regions',
                                                      'ratio'],
                                                  ascending=False)[
                                                    :filtering_threshold]
                                     .sort_values(by=['ratio',
                                                      'intersection_regions',
                                                      'pipeline_a_wt_regions',
                                                      'pipeline_b_regions'],
                                                  ascending=False)
                                     .reset_index(drop=True))
    rankings_intcount_filtered_df.index += 1

    if args.filter:
        peako_file_filtered = os.path.join(outdir,
                                           'peaKO-rankings-filtered.txt')
        rankings_intcount_filtered_df.to_csv(peako_file_filtered, index=True,
                                             sep='\t')
        logging.info('Wrote peaKO file filtered (top 50% by number of ' +
                     'intersecting regions).')

    # find and print motif's original vs. new ranking
    if jaspar_id or transcription_factor:
        if jaspar_id:
            logging.info('Using JASPAR ID ({}) to '.format(jaspar_id) +
                         'find ranks.')

            label = jaspar_id
            my_regex = "id" + jaspar_id
            column = "motif_ID"

        elif transcription_factor:
            logging.info('Using motif common name ({}) to find ranks.'.format(
                    transcription_factor))

            label = transcription_factor
            my_regex = "(?i)^" + transcription_factor + "$"  # case insensitive
            column = "alt_ID"

        # record motif's ranking as index (peaKO without filtering)
        peako_rank = (rankings_ratio_sorted_df[rankings_ratio_sorted_df[column]
                                               .str.contains(my_regex)]
                      .index[0])

        logging.info('PeaKO ranks {} in position {}.'
                     .format(label, peako_rank))

        try:
            peako_with_filtering_rank = (rankings_intcount_filtered_df[
                rankings_intcount_filtered_df[column]
                .str.contains(my_regex)]
                                         .index[0])
        except IndexError:
            peako_with_filtering_rank = 'NA (motif has been filtered out)'

        if args.filter:
            logging.info('PeaKO after 50% filtering by number of ' +
                         'intersections ranks {} in position {}.'
                         .format(label, peako_with_filtering_rank))
    else:
        logging.warning('User did not provide any motif IDs.')

    if extra:
        # print all rankings to file to use in plotting
        # the user can manually "cat" all output files afterwards
        logging.info('Generating CSV file to be used for plotting...')

        if jaspar_id:
            plotting_file = os.path.join(outdir,
                                         '{}-{}-centrimo-plot-data.txt'.format(
                                             transcription_factor, jaspar_id))
            # find DMC rank
            fe_df = read_and_sort_fisher_exact_centrimo_df(dmc_txt)
            dmc_rank = fe_df.index[fe_df["id"] == jaspar_id].tolist()[0]
            # find KOIN rank
            koin_rank = sorted_koin_df.index[
                sorted_koin_df["id"] == jaspar_id].tolist()[0]
            # find WT rank (from DMC sorted by E-value instead of
            # Fisher's E-value)
            wt_rank = sorted_dmc_df.index[
                sorted_dmc_df["id"] == jaspar_id].tolist()[0]
            # find number of peako motifs
            num_peako_motifs = len(rankings_ratio_sorted_df)
            # find number of peaks WT
            data_dmc_html = store_data_from_html(dmc_html)
            num_wt_peaks = data_dmc_html["sequence_db"]["count"]
            # find number of peaks KO
            num_ko_peaks = data_dmc_html["negative_sequence_db"]["count"]
            # find number of peaks KOIN
            data_koin_html = store_data_from_html(koin_html)
            num_koin_peaks = data_koin_html["sequence_db"]["count"]
            with open(plotting_file, 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter=',',
                                    lineterminator='\n')
                # Experiment, DMC, KOIN, WT, peaKO, peaKO_wo_filtering,
                # motifs, WT_peaks, KO_peaks, KOIN_peaks
                writer.writerow(["Experiment", "DMC", "KOIN", "WT", "peaKO",
                                 "peaKO_wo_filtering",
                                 "num_motifs_peaKO_wo_filt",
                                 "num_motifs_DMC", "num_motifs_KOIN",
                                 "WT_peaks", "KO_peaks", "KOIN_peaks"])
                writer.writerow([transcription_factor, dmc_rank, koin_rank,
                                 wt_rank, peako_with_filtering_rank,
                                 peako_rank, num_peako_motifs,
                                 number_of_motifs, number_of_koin_motifs,
                                 num_wt_peaks, num_ko_peaks,
                                 num_koin_peaks])

    logging.info('Finished! :)')


if __name__ == '__main__':
    main()


# TO DO:
# change output plotting file names to match requirements of plotting scripts
# when printing plotting file, use tuples to group column names and values

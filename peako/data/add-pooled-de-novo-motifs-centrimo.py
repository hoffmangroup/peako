#!/usr/bin/env python
from __future__ import print_function, division
from __future__ import absolute_import
import argparse
from bs4 import BeautifulSoup
import yaml
from collections import OrderedDict
import csv
import os
import re
import six
from six.moves import filter

parser = argparse.ArgumentParser(description='Add pooled MEME and DREME ' +
                                 'motif databases (Pipeline A, Pipeline B, ' +
                                 'and WT) to CentriMo calls. Must have ' +
                                 'already run pool-de-novo-motifs-before-' +
                                 'centrimo.py.')
parser.add_argument('-l', '--list', nargs='+',
                    help='list of CentriMo .html files', required=True)
parser.add_argument('-o', help='output directory for CentriMo', required=True)
args = parser.parse_args()


# PARAMETERS

centrimo_file_list = args.list
output_directory = args.o  # for CentriMo call output
# (new script will be created in $(pwd)/centrimo-neg_out)

# regexes to find dreme or meme motifs databases (.xml)
# and Pipeline A, Pipeline B, and WT files, and JASPAR db
regex = re.compile(r'dreme.xml|meme.xml')
koin_regex = re.compile(r'.*PipelineB.*')
dmc_regex = re.compile(r'.*PipelineA.*')
wt_regex = re.compile(r'.*WT.*')
db_regex = re.compile(r'.*JASPAR.*')


# FUNCTIONS

# Adapted from: http://stackoverflow.com/questions/13323976/
def parse_html(htmlfile):
    """
    Read an HTML file and return text between 'script' tags
    as string objects.
    """
    opened_file = open(htmlfile)
    html = opened_file.read()
    soup = BeautifulSoup(html, "html.parser")
    opened_file.close()
    script = soup.find('script')
    json_text = script.text.split('=', 1)[1].strip().rstrip(';')
    data = yaml.safe_load(json_text)  # loads into string objects
    return data


def change_to_pool_filename(filename):
    """Return corresponding full filepath of pooled motif file."""
    return (os.path.dirname(filename) + '/' +
            os.path.splitext(os.path.basename(filename))[0] +
            '-pool.xml')


# SCRIPT

cmd_dictionary = OrderedDict()
pooled_de_novo_motif_db_dictionary = OrderedDict()

for filename in centrimo_file_list:
    filebase = os.path.basename(filename)

    centrimo_data = parse_html(filename)
    # centrimo_data['cmd'][0] = "centrimo-neg"  # change from 'centrimo'
    centrimo_data['cmd'][6] = (output_directory + '/centrimo-' +
                               str(centrimo_file_list.index(filename)) + '/')

    # check filepath for regexes defined above
    # if present, then append meme-pool.xml and dreme-pool.xml to a
    # list or dictionary
    # XXX should simplify this section; lots of repeated code
    if dmc_regex.search(filename):
        cmd = [change_to_pool_filename(i) if regex.search(i)
               else i for i in centrimo_data['cmd']]
        pooled_de_novo_motif_db_dictionary['DMC'] = [
            change_to_pool_filename(i) for i in centrimo_data['cmd']
            if regex.search(i)]
        simple_file = 'DMC'
    elif koin_regex.search(filename):
        cmd = [change_to_pool_filename(i) if regex.search(i)
               else i for i in centrimo_data['cmd']]
        pooled_de_novo_motif_db_dictionary['KOIN'] = [
            change_to_pool_filename(i) for i in centrimo_data['cmd']
            if regex.search(i)]
        simple_file = 'KOIN'
    elif wt_regex.search(filename):
        cmd = [change_to_pool_filename(i) if regex.search(i)
               else i for i in centrimo_data['cmd']]
        pooled_de_novo_motif_db_dictionary['WT'] = [
            change_to_pool_filename(i) for i in centrimo_data['cmd']
            if regex.search(i)]
        simple_file = 'WT'
    else:
        cmd = centrimo_data['cmd']
        simple_file = (os.path.splitext(filebase)[0] +
                       str(centrimo_file_list.index(filename)))
        # remove extension

    cmd_dictionary[simple_file] = cmd


# insert motif databases (.xml files) from DMC, KOIN, WT runs
for key, cmd in six.iteritems(cmd_dictionary):
    if key == 'DMC':
        insert_motif_dbs = (pooled_de_novo_motif_db_dictionary['KOIN'] +
                            pooled_de_novo_motif_db_dictionary['WT'])
        for i in insert_motif_dbs:
            cmd.insert(cmd.index(list(filter(db_regex.match, cmd))[0]), i)
    elif key == 'KOIN':
        insert_motif_dbs = (pooled_de_novo_motif_db_dictionary['DMC'] +
                            pooled_de_novo_motif_db_dictionary['WT'])
        for i in insert_motif_dbs:
            cmd.insert(cmd.index(list(filter(db_regex.match, cmd))[0]), i)
    elif key == 'WT':
        insert_motif_dbs = (pooled_de_novo_motif_db_dictionary['DMC'] +
                            pooled_de_novo_motif_db_dictionary['KOIN'])
        for i in insert_motif_dbs:
            cmd.insert(cmd.index(list(filter(db_regex.match, cmd))[0]), i)
    else:
        insert_motif_dbs = (pooled_de_novo_motif_db_dictionary['KOIN'] +
                            pooled_de_novo_motif_db_dictionary['DMC'] +
                            pooled_de_novo_motif_db_dictionary['WT'])
        for i in insert_motif_dbs:
            cmd.insert(cmd.index(list(filter(db_regex.match, cmd))[0]), i)
    cmd_dictionary[key] = ' '.join(cmd)


# print commands to file
output_file = os.path.join('.', 'commands-pooled-de-novo.sh')

with open(output_file, 'w') as csv_file:
    writer = csv.writer(csv_file, delimiter='\t', lineterminator='\n\n')
    writer.writerow(["#!/usr/bin/env bash"])
    for _, cmd in cmd_dictionary.items():
        writer.writerow([cmd])

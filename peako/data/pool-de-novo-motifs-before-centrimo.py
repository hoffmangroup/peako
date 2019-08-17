#!/usr/bin/env python
from __future__ import print_function, division
from __future__ import absolute_import
import argparse
import os
from xml.etree import ElementTree as et
import logging
logging.basicConfig(format='%(asctime)s %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

parser = argparse.ArgumentParser(description='Modify ID of MEME and DREME ' +
                                 'de novo motifs to preserve uniqueness ' +
                                 'upon pooling DMC and KOIN pipelines.')
parser.add_argument('-d', help='full path of DMC directory')
parser.add_argument('-k', help='full path of KOIN directory')
parser.add_argument('-w', help='full path of WT directory')

args = parser.parse_args()

# dreme and meme files (and folders) can just be found based off
# of this because they always follow the same pattern
# (when output by centrimo)

meme_xml_dmc = args.d + '/meme_out/meme.xml'
dreme_xml_dmc = args.d + '/dreme_out/dreme.xml'

meme_xml_koin = args.k + '/meme_out/meme.xml'
dreme_xml_koin = args.k + '/dreme_out/dreme.xml'

meme_xml_wt = args.w + '/meme_out/meme.xml'
dreme_xml_wt = args.w + '/dreme_out/dreme.xml'

# PARAMETERS

output_directory_meme_dmc = os.path.dirname(meme_xml_dmc)
output_directory_meme_koin = os.path.dirname(meme_xml_koin)
output_directory_meme_wt = os.path.dirname(meme_xml_wt)

output_directory_dreme_dmc = os.path.dirname(dreme_xml_dmc)
output_directory_dreme_koin = os.path.dirname(dreme_xml_koin)
output_directory_dreme_wt = os.path.dirname(dreme_xml_wt)


# FUNCTIONS


def replace_meme_id(meme_root, suffix):
    for motif in meme_root.iter('motif'):
        new_id = motif.attrib['name'] + suffix
        motif.set('name', new_id)


def add_motifs_to_dct(dreme_root, label):
    for motif in dreme_root.iter('motif'):
        motif_id = motif.attrib['seq']
        if motif_id in dreme_dct:
            dreme_dct[motif_id].append(label)
        else:
            dreme_dct[motif_id] = [label]


def append_suffix_dreme_id(dreme_root):
    for motif in dreme_root[1].findall('motif'):
        # XXX assume that 'motifs' is second layer
        motif_id = motif.attrib['seq']
        if dreme_dct[motif_id][-1] == 'used':
            dreme_root[1].remove(motif)
        else:
            suffix = '-' + '-'.join(dreme_dct[motif_id])
            new_id = motif_id + suffix
            motif.set('seq', new_id)
            dreme_dct[motif_id].append('used')


# BUILD TREE
def get_root(xml):
    tree = et.parse(xml)
    return tree.getroot()


tree_meme_dmc = et.parse(meme_xml_dmc)
root_meme_dmc = tree_meme_dmc.getroot()

tree_meme_koin = et.parse(meme_xml_koin)
root_meme_koin = tree_meme_koin.getroot()

tree_meme_wt = et.parse(meme_xml_wt)
root_meme_wt = tree_meme_wt.getroot()

tree_dreme_dmc = et.parse(dreme_xml_dmc)
root_dreme_dmc = tree_dreme_dmc.getroot()

tree_dreme_koin = et.parse(dreme_xml_koin)
root_dreme_koin = tree_dreme_koin.getroot()

tree_dreme_wt = et.parse(dreme_xml_wt)
root_dreme_wt = tree_dreme_wt.getroot()

logging.info('Built trees and got roots!')

# Alter MEME motif labels
replace_meme_id(root_meme_dmc, '-DMC')
replace_meme_id(root_meme_koin, '-KOIN')
replace_meme_id(root_meme_wt, '-WT')
logging.info('Renamed MEME motifs.')

# Create dictionary of DREME motif labels
dreme_dct = {}
add_motifs_to_dct(root_dreme_dmc, 'DMC')
add_motifs_to_dct(root_dreme_koin, 'KOIN')
add_motifs_to_dct(root_dreme_wt, 'WT')

# Alter DREME motif labels
append_suffix_dreme_id(root_dreme_dmc)
append_suffix_dreme_id(root_dreme_koin)
append_suffix_dreme_id(root_dreme_wt)
logging.info('Renamed DREME motifs and removed duplicates.')

# Comment:
# XXX It would be nice to output everything to the same file, but then
# I would have to insert between motif tags. This is doable BUT I think
# it's too confusing because of the information under 'model', unless
# we remove 'model' and 'run_time' tags completely.

# could change 'last_mod_dat' attribute of first child (tag 'model') of root
tree_meme_dmc.write(output_directory_meme_dmc + '/meme-pool.xml')
tree_meme_koin.write(output_directory_meme_koin + '/meme-pool.xml')
tree_meme_wt.write(output_directory_meme_wt + '/meme-pool.xml')

tree_dreme_dmc.write(output_directory_dreme_dmc + '/dreme-pool.xml')
tree_dreme_koin.write(output_directory_dreme_koin + '/dreme-pool.xml')
tree_dreme_wt.write(output_directory_dreme_wt + '/dreme-pool.xml')
logging.info('Wrote new trees to files. Check output directories!')
